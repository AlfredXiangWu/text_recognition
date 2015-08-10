/****************************************************************************/
//
// Matlab interface file:  psmFit_gpu.cu
//
// Written 2/12 by N. Howe
//
// Fits part-structured model to binary image using approximate GPU algorithm.
// Written for CUDA architecture.
//
// [E,loc] = psmFit_gpu(model,bimg,[max_x max_y])
//
/****************************************************************************/

#include <math.h>
#include "mex.h"

#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree

#define LARGE 1e10
#define BLOCK 192
#define INFINITY 0x7f800000

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
// macro for computing minimum
#endif

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
// macro for computing maximum
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
// macro for computing squares
#endif

#ifndef BOUND
#define BOUND(v, l, h) (((v) > (l)) ? (((v) < (h)) ? (v) : (h)) : (l))
  /* macro for bounding a number by two others */
#endif

#ifndef ABS
#define ABS(a) (((a) < 0) ? -(a) : (a))
// macro for computing absolute value
#endif

#ifndef ROUND
//#define ROUND(a) ((int)(floor((a)+0.5)))
#define ROUND(a) (((a)<0) ? ((int)((a)-0.5)) : ((int)((a)+0.5)))
// macro for computing rounded value
#endif

#ifndef FLOOR
#define FLOOR(a) (((a)<0) ? (((int)(a))-1) : ((int)(a)))
// macro for computing floor value
#endif

#define errCheck(a,b) if (!(a)) mexErrMsgTxt((b));

texture<float, 2, cudaReadModeElementType> texRef;

/* Kernel to calculate the cumulative sum of the array on the GPU */
__global__ void 
gdt_row_noloc(float* out, int m, int n, int hzn, size_t p, float xm) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if ((i<m)&&(j<n)) {
    // we are in bounds
    float s = INFINITY;
    float t;
    for (int k = MAX(-j,-hzn); k <= MIN(n-j-1,hzn); k++) {
      t = tex2D(texRef,i+0.5f,(j+k)+0.5f)+k*k/xm;
      s = (s<t) ? s:t;
    }
    out[i+j*p] = s;
  }
}

/* Kernel to calculate the cumulative sum of the array on the GPU */
__global__ void 
gdt_col_noloc(float* out, int m, int n, int hzn, size_t p, float ym) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if ((i<m)&&(j<n)) {
    // we are in bounds
    float s = INFINITY;
    float t;
    for (int k = MAX(-i,-hzn); k <= MIN(m-i-1,hzn); k++) {
      t = tex2D(texRef,i+k+0.5f,j+0.5f)+k*k/ym;
      s = (s<t) ? s:t;
    }    
    out[i+j*p] = s;
  }
}


/* Kernel to calculate the cumulative sum of the array on the GPU */
__global__ void gdt_row(float* out, int* loc, 
           int m, int n, int hzn, size_t p, float xm) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if ((i<m)&&(j<n)) {
    // we are in bounds
    float s = tex2D(texRef,i+0.5f,j+0.5f);
    float t;
    int l = 0;
    for (int k = MAX(-j,-hzn); k <= MIN(n-j-1,hzn); k++) {
      t = tex2D(texRef,i+0.5f,(j+k)+0.5f)+k*k/xm;
      l = (s<=t) ? l:k;
      s = (s<=t) ? s:t;
    }
    out[i+j*p] = s;
    loc[i+j*m] = j+l;
  }
}


/* Kernel to calculate the cumulative sum of the array on the GPU */
__global__ void gdt_col(float* out, int* loc, 
           int m, int n, int hzn, size_t p, float ym) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if ((i<m)&&(j<n)) {
    // we are in bounds
    float s = tex2D(texRef,i+0.5f,j+0.5f);
    float t;
    int l = 0;
    for (int k = MAX(-i,-hzn); k <= MIN(m-i-1,hzn); k++) {
      t = tex2D(texRef,i+k+0.5f,j+0.5f)+k*k/ym;
      l = (s<=t) ? l:k;
      s = (s<=t) ? s:t;
    }    
    out[i+j*p] = s;
    loc[i+j*m] = i+l;
  }
}


/* Kernel to compute the final location indices */
__global__ void loc_replace(int* loc, int *loc1, int* loc2, int m, int n) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if ((i < m)&&(j < n)) {
        // Don't continue if past end of array  
        if (loc2[i+j*m] < m) {
            loc[i+j*m] = loc1[loc2[i+j*m]+m*j];
            //loc[i+j*m] = loc1[i+j*m];
        } else {
            loc[i+j*m] = -2;
        }
    }
}

/* Kernel adds two buffers: c = a+b */
__global__ void add_kernel(float* a, float *b, float *c, int m, int n, int p) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if ((i < m)&&(j < n)) {
        // Don't continue if past end of array  
        c[i+j*p] = a[i+j*p]+b[i+j*p];
    }
}

/* Kernel multiplies two buffers: c = a*b */
__global__ void mult_kernel(float* a, float *b, float *c, int m, int n, int p) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if ((i < m)&&(j < n)) {
        // Don't continue if past end of array  
        c[i+j*p] = a[i+j*p]*b[i+j*p];
    }
}

/* Kernel copies values b <- a*/
__global__ void copy_kernel(float* a, float *b, int m, int n, int p) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if ((i < m)&&(j < n)) {
        // Don't continue if past end of array  
        b[i+j*p] = a[i+j*p];
    }
}

/* Kernel assigns uniform value */
__global__ void assign_kernel(float* a, float v, int m, int n, int p) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if ((i < m)&&(j < n)) {
        // Don't continue if past end of array  
        a[i+j*p] = v;
    }
}

/* Handle any CUDA errors */
void checkForCudaError(cudaError_t err) {
    if (err != cudaSuccess) {
        /* print message about this error */
        mexPrintf("CUDA error:  ");
        mexErrMsgTxt(cudaGetErrorString(err));
    }
}

/* Adds two buffers on GPU:  c = a+b */
void gpu_add(float* a, float *b, float *c, int m, int n, int p) {
    dim3 dimGrid((m+31)/32,n);
    dim3 dimBlock(32); 
    add_kernel<<<dimGrid,dimBlock>>>(a, b, c, m, n, p);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    checkForCudaError(err);
}

/* Multiplies two buffers on GPU:  c = a*b */
void gpu_mult(float* a, float *b, float *c, int m, int n, int p) {
    dim3 dimGrid((m+31)/32,n);
    dim3 dimBlock(32); 
    mult_kernel<<<dimGrid,dimBlock>>>(a, b, c, m, n, p);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    checkForCudaError(err);
}

/* Assigns to a buffer on GPU */
void gpu_add(float* a, float v, int m, int n, int p) {
    dim3 dimGrid((m+31)/32,n);
    dim3 dimBlock(32); 
    assign_kernel<<<dimGrid,dimBlock>>>(a, v, m, n, p);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    checkForCudaError(err);
}

/* Kernel helps perform translation */
__global__ void tran_kernel(float* dt1, float *dt2, int m, int n, int p, 
        float maxd, int x1, int y1, float s0, float s1, float s2, float s3,
        float t0, float t1, float t2, float t3) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    // Don't continue if past end of array  
    if ((i < m)&&(j < n)) {
        if ((j+x1 < 1)||(i+y1<1)||(j+x1+2>=n)||(i+y1+2>=m)) {
            dt2[i+j*p] = maxd;
        } else {
        dt2[i+j*p] = t0*s0*dt1[i+y1-1+p*(j+x1-1)]
            +t0*s1*dt1[i+y1-1+p*(j+x1)]
            +t0*s2*dt1[i+y1-1+p*(j+x1+1)]
            +t0*s3*dt1[i+y1-1+p*(j+x1+2)]
            +t1*s0*dt1[i+y1+p*(j+x1-1)]
            +t1*s1*dt1[i+y1+p*(j+x1)]
            +t1*s2*dt1[i+y1+p*(j+x1+1)]
            +t1*s3*dt1[i+y1+p*(j+x1+2)]
            +t2*s0*dt1[i+y1+1+p*(j+x1-1)]
            +t2*s1*dt1[i+y1+1+p*(j+x1)]
            +t2*s2*dt1[i+y1+1+p*(j+x1+1)]
            +t2*s3*dt1[i+y1+1+p*(j+x1+2)]
            +t3*s0*dt1[i+y1+2+p*(j+x1-1)]
            +t3*s1*dt1[i+y1+2+p*(j+x1)]
            +t3*s2*dt1[i+y1+2+p*(j+x1+1)]
            +t3*s3*dt1[i+y1+2+p*(j+x1+2)];
        }
    }
}

/****************************************************************************/
//
// helper function to do matrix translation --
// uses interpolation for subpixel accuracy
//
void tranDT(float *dt1, float *dt2, int m, int n, int p, 
     float x, float y, float maxd) {
  int x1 = FLOOR(x);
  int y1 = FLOOR(y);
  float xf = x-x1;
  float yf = y-y1;
  float t0 = (((2-yf)*yf-1)*yf)/2;
  float t1 = ((3*yf-5)*yf*yf+2)/2;
  float t2 = (((4-3*yf)*yf+1)*yf)/2;
  float t3 = ((yf-1)*yf*yf)/2;
  //mexPrintf("%f %f %f %f\n",t0,t1,t2,t3);
  float s0 = (((2-xf)*xf-1)*xf)/2;
  float s1 = ((3*xf-5)*xf*xf+2)/2;
  float s2 = (((4-3*xf)*xf+1)*xf)/2;
  float s3 = ((xf-1)*xf*xf)/2;
  //mexPrintf("%f %f %f %f\n",s0,s1,s2,s3);

    dim3 dimGrid((m+31)/32,n);
    dim3 dimBlock(32); 
    //mexPrintf("In tranDT; p = %d\n",p);
    tran_kernel<<<dimGrid,dimBlock>>>(dt1, dt2, m, n, p/sizeof(float),
        maxd, x1, y1, s0, s1, s2, s3, t0, t1, t2, t3);
    //copy_kernel<<<dimGrid,dimBlock>>>(dt1, dt2, m, n);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    checkForCudaError(err);
    //mexPrintf("tran_kernel finished\n");
}

/****************************************************************************/
//
// helper function to do 2d minimum convolution
//
void psdt2d(float *cost, float *out, int *loc, int m, int n, int p,
            int xhzn, int yhzn, float xm, float ym, cudaDeviceProp pdev,
            float *scratch, float *scratch2, int *lscratch) {
    int *loc2 = loc+(m*n);
    int pf = p/sizeof(float);
    //mexPrintf("Pitch:  %d -> %d\n",p,pf);

    /* Set up execution configuration based on m rows */
    //int m0 = m%pdev.maxThreadsDim[0]+1;
    dim3 dimGrid1((m+BLOCK-1)/BLOCK,n);
    dim3 dimBlock1(BLOCK); 
    // m threads for row-wise calculation (m rows)

    /* copy original cost to scratch array */
    cudaError_t err = cudaMemcpy2D(scratch2,p,cost,p,m*sizeof(float),n,
        cudaMemcpyDeviceToDevice);
    checkForCudaError(err);
    //mexPrintf("Scratch set up.\n");
            
    /* Call row function on GPU */
    cudaBindTexture2D(0,texRef,scratch2,m,n,p);
    err = cudaGetLastError();
    checkForCudaError(err);
    //mexCallMATLAB(0,0,0,0,"tic");
    if (loc) {
        gdt_row<<<dimGrid1,dimBlock1>>>(scratch, lscratch, m, n, xhzn, pf, xm);
    } else {
        gdt_row_noloc<<<dimGrid1,dimBlock1>>>(scratch, m, n, xhzn, pf, xm);
    }
    cudaThreadSynchronize();
    //mexCallMATLAB(0,0,0,0,"toc");
    err = cudaGetLastError();
    checkForCudaError(err);
    cudaUnbindTexture(texRef);  
    err = cudaGetLastError();
    checkForCudaError(err);
    //mexPrintf("Row computation complete.\n");

    /* Set up execution configuration based on n columns */
    //int n0 = n%pdev.maxThreadsDim[0]+1;
    dim3 dimGrid2((m+BLOCK-1)/BLOCK,n);
    dim3 dimBlock2(BLOCK); 
    // n threads for column-wise calculation (n columns)
     
    /* Call column function on GPU */
    cudaBindTexture2D(0,texRef,scratch,m,n,p);
    err = cudaGetLastError();
    checkForCudaError(err);
    //mexCallMATLAB(0,0,0,0,"tic");
    if (loc) {
        gdt_col<<<dimGrid2,dimBlock2>>>(out, loc2, m, n, yhzn, pf, ym);
    } else {
        gdt_col_noloc<<<dimGrid2,dimBlock2>>>(out, m, n, yhzn, pf, ym);
    }
    cudaThreadSynchronize();
    //mexCallMATLAB(0,0,0,0,"toc");
    err = cudaGetLastError();
    checkForCudaError(err);
    cudaUnbindTexture(texRef);  
    err = cudaGetLastError();
    checkForCudaError(err);
    //mexPrintf("Column computation complete.\n");

    if (loc) {
        /* Set up execution configuration based on m x n pixels */     
        dim3 dimGrid3((m+31)/32,n);
        dim3 dimBlock3(32); 
        //mexPrintf("Blocks:  %d %d %d\n",dimBlock3.x,dimBlock3.y,dimBlock3.z);
        //mexPrintf("Grid:  %d %d %d\n",dimGrid3.x,dimGrid3.y,dimGrid3.z);

        /* Call loc function on GPU */
        loc_replace<<<dimGrid3,dimBlock3>>>(loc, lscratch, loc2, m, n);
        cudaThreadSynchronize();
        err = cudaGetLastError();
        checkForCudaError(err);
        //mexPrintf("Replacement done.\n");
    }
}

/****************************************************************************/
//
// gateway driver to call the distance transform code from matlab
//
// This is the matlab entry point
void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i, j, c, nrow, ncol, nmpt, xhzn, yhzn;
  int mxfn, myfn, mrfn, mvxfn, mvyfn, mpfn, mchfn;
  int root, curt, prnt, ch, dtarrbase, dtarrlim;
  int *mp, *mnch, *mdone, *mndesc, *dtarrhist;
  int *loc, *loc2, *lscratch, **locarr, **locarr2;
  char *pbimgc;
  double maxd;
  double *out, *scl;
  double *mx, *my, *mr, *mvx, *mvy, *pbimg;
  double **mch, **locp; 
  float *scratch, *scratch2, *dt1, *tranarr;
  float **dtarr, **dtarrp;
  float *imgb;
  mxArray *cell, *field;
  size_t p;
  cudaError_t err;

  //mexCallMATLAB(0,0,0,0,"tic");

  // check for proper number and size of arguments
  errCheck(nrhs == 3,"Arguments:  model, bimg, hzn");
  //errCheck(nlhs <= 3,"Outputs:  dt, [loc], [scale]");
  //errCheck(nlhs <= 1,"Outputs:  dt (loc disabled/not implemented)");
  errCheck(mxIsStruct(prhs[0]),"model must be struct.");
  mxfn = mxGetFieldNumber(prhs[0],"x");
  errCheck(mxfn>=0,"Model must have field named x.");
  myfn = mxGetFieldNumber(prhs[0],"y");
  errCheck(myfn>=0,"Model must have field named y.");
  mrfn = mxGetFieldNumber(prhs[0],"r");
  errCheck(mrfn>=0,"Model must have field named r.");
  mvxfn = mxGetFieldNumber(prhs[0],"vx");
  errCheck(mvxfn>=0,"Model must have field named vx.");
  mvyfn = mxGetFieldNumber(prhs[0],"vy");
  errCheck(mvyfn>=0,"Model must have field named vy.");
  mpfn = mxGetFieldNumber(prhs[0],"parent");
  errCheck(mpfn>=0,"Model must have field named parent.");
  mchfn = mxGetFieldNumber(prhs[0],"children");
  errCheck(mchfn>=0,"Model must have field named children.");
  errCheck(mxIsDouble(prhs[1])||mxIsLogical(prhs[1])||mxIsUint8(prhs[1]),
           "Second argument must be binary image type.");
  errCheck(mxGetNumberOfDimensions(prhs[1])==2,
           "Second argument must be binary image type.");
  errCheck(mxIsDouble(prhs[2])&&!mxIsComplex(prhs[2]),
             "hzn must be real double.");
  if (mxGetNumberOfElements(prhs[2])==1) {
    mexErrMsgTxt("Must now specify [x y] horizon in pixels.");
    ////xhzn = yhzn = (int)(*mxGetPr(prhs[2]));
    //maxd = *mxGetPr(prhs[2]);
    //xhzn = yhzn = CEIL(sqrt(maxd/2));
  } else if (mxGetNumberOfElements(prhs[2])==2) {
    xhzn = (int)(*mxGetPr(prhs[2]));
    yhzn = (int)(*(mxGetPr(prhs[2])+1));
    maxd = xhzn*xhzn+yhzn*yhzn;
  } else {
    mexErrMsgTxt("hzn must be scalar or 2-vector.");
  }
  //mexPrintf("Horizons:  %d & %d\n",xhzn,yhzn);

  //mexPrintf("maxd = %g\n",maxd);
  nmpt = (int)mxGetNumberOfElements(prhs[0]);
  nrow = (int)mxGetM(prhs[1]);
  ncol = (int)mxGetN(prhs[1]);
  //mexPrintf("Inputs checked.\n");

  // transcribe model variables
  root = -1;
  mx = (double*)mxMalloc(nmpt*sizeof(double));
  my = (double*)mxMalloc(nmpt*sizeof(double));
  mr = (double*)mxMalloc(nmpt*sizeof(double));
  mvx = (double*)mxMalloc(nmpt*sizeof(double));
  mvy = (double*)mxMalloc(nmpt*sizeof(double));
  mp = (int*)mxMalloc(nmpt*sizeof(int));
  mch = (double**)mxMalloc(nmpt*sizeof(double*));
  mnch = (int*)mxMalloc(nmpt*sizeof(int));
  for (int i = 0; i < nmpt; i++) {
    mx[i] = *mxGetPr(mxGetFieldByNumber(prhs[0],i,mxfn));
    errCheck(ABS(ROUND(mx[i])) < ncol,"Model displacement larger than image.");
    my[i] = *mxGetPr(mxGetFieldByNumber(prhs[0],i,myfn));
    errCheck(ABS(ROUND(my[i])) < nrow,"Model displacement larger than image.");
    mr[i] = *mxGetPr(mxGetFieldByNumber(prhs[0],i,mrfn));
    mvx[i] = *mxGetPr(mxGetFieldByNumber(prhs[0],i,mvxfn));
    mvy[i] = *mxGetPr(mxGetFieldByNumber(prhs[0],i,mvyfn));
    mp[i] = (int)(*mxGetPr(mxGetFieldByNumber(prhs[0],i,mpfn)));
    field = mxGetFieldByNumber(prhs[0],i,mchfn);
    mch[i] = mxGetPr(field);
    mnch[i] = (int)(mxGetNumberOfElements(field));
    if (mp[i]==0) {
      errCheck(root<0,"Model must have a single root.");
      root = i;
    }
  }
  errCheck(root>=0,"Model must have a root.");
  //mexPrintf("Model checked.\n");

  //for (int i = 0; i < nmpt; i++) {
  //mexPrintf("%d: (%f,%f)\n",i,mvx[i],mvy[i]);
  //}

  // allocate output space
  plhs[0] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  out = mxGetPr(plhs[0]);
  if (nlhs>=2) {
    locp = (double**)mxMalloc(nrow*ncol*sizeof(double*));
    plhs[1] = mxCreateCellMatrix(nrow,ncol);
    for (i = 0; i < nrow*ncol; i++) {
      cell = mxCreateDoubleMatrix(2,nmpt,mxREAL);
      mxSetCell(plhs[1],i,cell);
      locp[i] = mxGetPr(cell);
    }
    loc = (int*)mxMalloc(2*nrow*ncol*sizeof(int));
    loc2 = loc+nrow*ncol;
    if (nlhs>=3) {
      plhs[2] = mxCreateDoubleMatrix(nrow,ncol,mxREAL);
      scl = mxGetPr(plhs[2]);
    }
  } else {
    loc = 0;
    //loc2 = 0;
    locp = 0;
  }
  //mexPrintf("Output allocated.\n");

  // get gpu device information
  int  ndev;
  cudaGetDeviceCount(&ndev);
  cudaThreadSynchronize();
  //mexPrintf("There are %d GPUs.\n",ndev);
  if (ndev < 1) {
    mexErrMsgTxt("No CUDA device found.");
  }
  cudaDeviceProp pdev;
  cudaGetDeviceProperties(&pdev,0);
  cudaThreadSynchronize();
  //mexPrintf("GPU device found.\n");
  err = cudaDeviceReset();
  checkForCudaError(err);

  // allocate GPU space
  err = cudaMallocPitch((void **)&dt1,&p,sizeof(float)*nrow,ncol);
  //err = cudaMalloc((void **)&dt1,sizeof(float)*nrow*ncol);
  checkForCudaError(err);
  err = cudaMallocPitch((void **)&scratch,&p,sizeof(float)*nrow,ncol);
  //err = cudaMalloc((void **)&scratch,sizeof(float)*nrow*ncol);
  checkForCudaError(err);
  err = cudaMallocPitch((void **)&scratch2,&p,sizeof(float)*nrow,ncol);
  //err = cudaMalloc((void **)&scratch2,sizeof(float)*nrow*ncol);
  checkForCudaError(err);
  if (loc) {
    err = cudaMalloc((void **)&lscratch,sizeof(int)*nrow*ncol);
    checkForCudaError(err);
  } else {
    lscratch = 0;
  }
  /*err = cudaMalloc( (void **) &v,sizeof(int)*nrow*ncol);   // Allocate memory
  checkForCudaError(err);
  err = cudaMalloc( (void **) &z,sizeof(float)*(nrow+1)*(ncol+1));
  checkForCudaError(err);*/
  err = cudaMallocPitch( (void **) &tranarr, &p, sizeof(float)*nrow, ncol); 
  // stores translation
  checkForCudaError(err);
  //mexPrintf("GPU storage allocated.\n");
       
  // set up textures
  texRef.normalized = false;
  texRef.filterMode = cudaFilterModeLinear;
  texRef.addressMode[0] = cudaAddressModeWrap;
  texRef.addressMode[1] = cudaAddressModeWrap;
       
  // transcribe image into dt1
  pbimg = mxGetPr(prhs[1]);
  pbimgc = (char*)pbimg;
  imgb = (float*)mxMalloc(nrow*ncol*sizeof(float));
  if (mxIsDouble(prhs[1])) {
    for (int i = 0; i <nrow*ncol; i++) {
      if (pbimg[i]) {
        imgb[i] = 0;
      } else {
        imgb[i] = (float)maxd;
      }
    }
  } else {
    for (int i = 0; i <nrow*ncol; i++) {
      if (pbimgc[i]) {
        imgb[i] = 0;
      } else {
        imgb[i] = (float)maxd;
      }
    }
  }

  err = cudaMemcpy2D(dt1,p,imgb,nrow*sizeof(float),nrow*sizeof(float),ncol,
    cudaMemcpyHostToDevice);
  checkForCudaError(err);
  //mexPrintf("Image transcribed.\n");

  // allocate other variables
  mdone = (int*)mxMalloc(nmpt*sizeof(int));
  mndesc = (int*)mxMalloc(nmpt*sizeof(int));
  dtarr = (float**)mxMalloc(nmpt*sizeof(float*));  // stores intermediate dt
  dtarrp = (float**)mxMalloc(nmpt*sizeof(float*));  // pointer to storage
  dtarrhist = (int*)mxMalloc(nmpt*sizeof(int));  // stores historical size
  locarr = (int**)mxMalloc(2*nmpt*sizeof(int*));  // stores intermediate loc
  locarr2 = locarr+nmpt;
  for (i = 0; i < nmpt; i++) {
    mdone[i] = 0;
    mndesc[i] = 1;
    dtarrp[i] = 0;
  }
  if (loc) {
    for (i = 0; i < nmpt; i++) {
      err = cudaMalloc((void **)locarr+i,sizeof(int)*2*nrow*ncol);
      checkForCudaError(err);
      locarr2[i] = locarr[i]+nrow*ncol;
    }
  } else {
    for (i = 0; i < nmpt; i++) {
      locarr[i] = 0;
      locarr2[i] = 0;
    }
  }
  //mexPrintf("Other variables allocated.\n");

  // compute distance transform on input image
  psdt2d(dt1,dt1,0,nrow,ncol,p,xhzn,yhzn,1,1,pdev,scratch,scratch2,lscratch);
  //mexPrintf("Initial distance transform.\n");

  /*
  // debug test
  //tranDT(dt1, scratch, nrow, ncol, p, 5.0, 0.0, maxd);
  //tranDT(scratch, dt1, nrow, ncol, p, 5.0, 0.0, maxd);

  // DEBUG:  view dt image
  //mexPrintf("Pitches:  %d, %d\n",sizeof(float)*nrow,p);
  err = cudaMemcpy2D(imgb,sizeof(float)*nrow,dt1,p,sizeof(float)*nrow,ncol, 
    cudaMemcpyDeviceToHost);
  checkForCudaError(err);
  for (i = 0; i < nrow*ncol; i++) {
    //mexPrintf("%d:  %f\n",i,imgb[i]);
    out[i] = double(imgb[i]);
  }
  //mexPrintf("Debug stop.\n");
  return;*/

  // run computation
  curt = root;
  dtarrlim = 0;
  dtarrbase = 0;
  dtarrp[root] = dt1;
  //mexPrintf("Start:  node %d, done %d/%d\n",curt,mdone[curt],mnch[curt]);
  //mexCallMATLAB(0,0,0,0,"tic");
  //mexCallMATLAB(0,0,0,0,"toc");
  while (curt != -1) {
    //mexPrintf("Looping:  node %d, done %d/%d\n",curt,mdone[curt],mnch[curt]);

    // set base point
    if (mdone[curt]==0) {
      dtarrhist[curt] = dtarrbase;
    }

    // figure out if we need to descend to children
    if (mdone[curt]<mnch[curt]) {
      // figure out current child
      ch = (int)(mch[curt][mdone[curt]])-1;  // Matlab indices start from 1
        
      //if (dtarrp[ch]==0) {
      //mexPrintf("Set up memory to do child %d out of %d for %d (in %d)\n",
      //          mdone[curt],mnch[curt],curt,dtarrbase);
      
      // we haven't yet computed the current child; set up to do that now
      ch = (int)(mch[curt][mdone[curt]])-1;  // Matlab indices start from 1
      
      // allocate new dt storage space if necessary
      if (dtarrlim==dtarrbase) {
        err = cudaMallocPitch((void **)dtarr+dtarrlim,&p,
            sizeof(float)*nrow,ncol);
        checkForCudaError(err);
        dtarrlim++;
        //mexPrintf("New dt allocation.\n");
      }      

      // copy dt1
      err = cudaMemcpy2D(dtarr[dtarrbase],p,dt1,p,nrow*sizeof(float),ncol, 
          cudaMemcpyDeviceToDevice);
      checkForCudaError(err);
      //mexPrintf("Copied dt1.\n");
                 
      // set up storage for this child
      dtarrp[ch] = dtarr[dtarrbase];
      dtarrbase++;
      
      // move on, but when we return, remember that this child is done
      mdone[curt]++;
      curt = ch;
      //} else {
      //mexPrintf("Advance to next child for %d\n",curt);
      //
      // this child is complete; advance to the next if there is one
      //mdone[curt]++;
      //}
    } else {
      //mexPrintf("Assemble result for %d (%d children)\n",curt,mnch[curt]);

      // all children complete; assemble result
      for (c = 0; c < mnch[curt]; c++) {
        ch = (int)(mch[curt][c])-1;  // Matlab indices start from 1
        //mexPrintf("Beginning translation.\n");
        //mexPrintf("  %d translate (%d,%d), var (%f,%f)\n",
        //          ch,ROUND(mx[ch]),ROUND(my[ch]),mvx[ch],mvy[ch]);
        //mexPrintf("  %d actual (%g,%g)\n",ch,mx[ch],my[ch]);
        errCheck(dtarrp[ch]!=0,"Null pointer");
        //mexPrintf("Before:  %f; %f\n",dtarrp[ch][0],dtarrp[ch][1]);
        mxAssert((ch>=0)&&(ch<nmpt),"Out of bounds");
        //mexPrintf("mndesc[%d] = %d\n",ch,mndesc[ch]);
        //mexPrintf("Translation!\n");
        //mexCallMATLAB(0,0,0,0,"tic");
        tranDT(dtarrp[ch],tranarr,nrow,ncol,p,mx[ch],my[ch],maxd*mndesc[ch]);
        //mexPrintf("Translation complete.\n");
        //mexCallMATLAB(0,0,0,0,"toc");
        //mexPrintf("After:  %f; %f\n",dtarrp[ch][0],dtarrp[ch][1]);
        //mexPrintf("Beginning distance transform.\n");
        psdt2d(tranarr,dtarrp[ch],locarr[ch],nrow,ncol,p,xhzn,yhzn,
          mvx[ch],mvy[ch],pdev,scratch,scratch2,lscratch);
        //mexCallMATLAB(0,0,0,0,"tic");
        //mexPrintf("Beginning addition.\n");
        gpu_add(dtarrp[curt],dtarrp[ch],dtarrp[curt],nrow,ncol,p/sizeof(float));
        //mexCallMATLAB(0,0,0,0,"toc");
      }

      // now go back up to the parent
      //mexPrintf("  %d %d\n",curt,mp[curt]-1);
      if (mp[curt]>0) {
        mxAssert((curt>=0)&&(curt<nmpt),"Out of bounds");
        mxAssert((mp[curt]-1>=0)&&(mp[curt]-1<nmpt),"Out of bounds");
        mndesc[mp[curt]-1] += mndesc[curt];
      }
      dtarrbase = dtarrhist[curt];
      curt = mp[curt]-1;  // Matlab indices start from 1
    }
  }
  //mexPrintf("Main computation complete.\n");
  //mexCallMATLAB(0,0,0,0,"toc");
  //mexCallMATLAB(0,0,0,0,"tic");

  // do final translation and distance transform, and copy result back
  tranDT(dt1,tranarr,nrow,ncol,p,mx[root],my[root],maxd*nmpt);
  psdt2d(tranarr,dt1,locarr[root],nrow,ncol,p,xhzn,yhzn,mvx[root],mvy[root],
         pdev,scratch,scratch2,lscratch);
  //mexPrintf("Copying to host.\n");
  err = cudaMemcpy2D(imgb,nrow*sizeof(float),dt1,p,nrow*sizeof(float),ncol, 
      cudaMemcpyDeviceToHost);
  checkForCudaError(err);
  for (i = 0; i < nrow*ncol; i++) {
    //mexPrintf("%d:  %f\n",i,imgb[i]);
    out[i] = double(imgb[i]);
  }
  //mexPrintf("Final translation complete.\n");

  // assemble loc information if required
  prnt = -1;
  if (loc) {
    for (i = 0; i < nmpt; i++) {
      mdone[i] = 0;
    }
    curt = root;

    while (curt != -1) {
      // first time through:
      if (mdone[curt] == 0) {
        int y = ROUND(my[curt]);
        int x = ROUND(mx[curt]);
        //mexPrintf("Node %d offset (%d,%d)\n",curt+1,x,y);
        
        // copy results for this node from gpu to cpu
        err = cudaMemcpy(loc, locarr[curt], 2*sizeof(int)*nrow*ncol, 
            cudaMemcpyDeviceToHost);
        checkForCudaError(err);
        
        // now copy to matlab result cells
        for (i = 0; i < nrow; i++) {
          for (j = 0; j < ncol; j++) {
            if (prnt == -1) {
              locp[i+j*nrow][2*curt] = BOUND(j+x,0,ncol-1);
              locp[i+j*nrow][2*curt+1] = BOUND(i+y,0,nrow-1);
              //if ((i==30)&&(j==71)) 
              //mexPrintf("(31,72): %d,%d + %d,%d = %f,%f\n",i+1,j+1,x,y,
              //          locp[i+j*nrow][2*curt]+1,locp[i+j*nrow][2*curt+1]+1);
            } else {
              int jj = (int)BOUND(locp[i+j*nrow][2*prnt],0,ncol-1);
              int ii = (int)BOUND(locp[i+j*nrow][2*prnt+1],0,nrow-1);
              locp[i+j*nrow][2*curt] = loc[ii+nrow*jj]+x;
              locp[i+j*nrow][2*curt+1] = loc2[ii+nrow*jj]+y;
              //if ((i==30)&&(j==71)) 
              //mexPrintf("(31,72): %d,%d ->  %d,%d + %d,%d = %f,%f\n",
              //          jj+1,ii+1,loc[ii+nrow*jj]+1,
              //          loc2[ii+nrow*jj]+1,
              //          x,y,locp[i+j*nrow][2*curt]+1,
              //          locp[i+j*nrow][2*curt+1]+1);
            }
          }
        }
      }

      //if (mdone[curt] == 0) {
      //  // copy results for this node
      //for (i = 0; i < nrow; i++) {
      //for (j = 0; j < ncol; j++) {
      //    //mexPrintf("%d %d:%d+%f\n",i,2*curt,locarr[curt][loc[i]],mx[curt]);
      //    mxAssert((2*curt>=0)&&(2*curt+1)<2*nmpt,"Bounds");
      //  locp[i+j*nrow][2*curt] = BOUND(locarr[curt][loc[i]]+x,0,ncol-1)+1;  
      //    // +1 for Matlab
      //locp[i+j*nrow][2*curt+1] = BOUND(locarr2[curt][loc[i]]+y,0,nrow-1+1);   
      //  // +1 for Matlab
      //loc[i+j*nrow] = (int)(locp[i][2*curt+1]+nrow*locp[i][2*curt]-nrow-1);
      //  // undo +1 for Matlab
      //    if ((i==30*nrow+71)&&(curt==13)) 
      //    mexPrintf("(30,71): %f %g\n",locp[i][2*curt],locp[i][2*curt+1]);
      //  }
      //}
      //}

      // descend to children if any left
      if (mdone[curt]<mnch[curt]) {
        // reset current level for second children
        //if (mdone > 0) {
        //for (i = 0; i < nrow*ncol; i++) {
        //  loc[i] = (int)(locp[i][2*curt-1]+nrow*locp[i][2*curt-2]);
        //}
        //}

        // descend to child
        prnt = curt;
        ch = (int)(mch[curt][mdone[curt]])-1;  // Matlab indices start from 1
        mdone[curt]++;
        curt = ch;
      } else {
        // come back up
        curt = prnt;  // Matlab indices start from 1
        prnt = mp[curt]-1;
      }
    }

    // add +1 to all indices for Matlab
    for (i = 0; i < nrow*ncol; i++) {
      for (j = 0; j < 2*nmpt; j++) {
        locp[i][j]++;
      }
    }
    //mexPrintf("loc assembled\n");
  }

  // compute scale if needed
  if (nlhs>=3) {
    //mexPrintf("Computing scale.\n");
    for (i = 0; i < nrow*ncol; i++) {
      double tmp = 0;
      
      for (j = 0; j < nmpt; j++) {
        if (mp[j]>0) {
          tmp += sqrt(SQR(locp[i][2*j]-locp[i][2*mp[j]-2])+
                      SQR(locp[i][2*j+1]-locp[i][2*mp[j]-1]));
        }
      }
      scl[i] = tmp;
    }
  }

  // free stuff
  for (i = 0; i < dtarrlim; i++) {
    cudaFree(dtarr[i]);
  }
  if (loc) {
    for (i = 0; i < nmpt; i++) {
      cudaFree(locarr[i]);
    }
  }
  mxFree(imgb);
  mxFree(locarr);
  mxFree(dtarrhist);
  mxFree(dtarrp);
  mxFree(dtarr);
  mxFree(mndesc);
  mxFree(mdone);
  cudaFree(tranarr);
  //cudaFree(z);
  //cudaFree(v);
  if (locp) {
    cudaFree(lscratch);
  }
  cudaFree(scratch2);
  cudaFree(scratch);
  cudaFree(dt1);
  if (nlhs==2) {
    mxFree(loc);
    mxFree(locp);
  }
  mxFree(mnch);
  mxFree(mch);
  mxFree(mp);
  mxFree(mvy);
  mxFree(mvx);
  mxFree(mr);
  mxFree(my);
  mxFree(mx);
  //mexPrintf("Freed.\n");  
  //mexCallMATLAB(0,0,0,0,"toc");
}
