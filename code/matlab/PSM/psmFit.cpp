/****************************************************************************/
//
// Matlab interface file:  psmFit.cpp
//
// Written 2/12 by N. Howe.
// Fits a part-structured model to a binary image.
// Uses bicubic interpolation for subpixel accuracy on translation.
//
// [dtsq,loc] = psmFit(model,bimg,maxd)
//
/****************************************************************************/

#include <math.h>
#include "mex.h"
#include "mymex.h"
#include "macros.h"

#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree

#define LARGE 1e10

/****************************************************************************/
//
// helper function to do 1d minimum convolution
//
void gdt1d(double *cost, double *out, int *loc, int *v, double *z, int sz_x) {
  int k, q;
  double s;

  // set up
  k = 0;
  v[0] = 0;
  z[0] = -mxGetInf();
  z[1] = mxGetInf();

  // compute
  for (q=1; q<sz_x; q++) {
    mxAssert(q<sz_x,"Q out of range");
    mxAssert(k<sz_x,"v out of range");
    mxAssert(v[k]+1<sz_x,"cost out of range");
    s = ((cost[q]+q*q)-(cost[v[k]]+v[k]*v[k]))/(2*(q-v[k]));  // intercept
    while (s<=z[k]) {
      k = k-1;
      mxAssert(q<sz_x,"Q out of range");
      mxAssert(k<sz_x,"v out of range");
      mxAssert(v[k]+1<sz_x,"cost out of range");
      s = ((cost[q]+q*q)-(cost[v[k]]+v[k]*v[k]))/(2*(q-v[k]));
    }
    k = k+1;
    mxAssert(k<sz_x,"k out of range");
    v[k] = q;
    z[k] = s;
    z[k+1] = mxGetInf();
  }
  k = 0;
  for (q=0; q<sz_x; q++) {
    while (z[k+1]<q) {
      k++;
    }
    if (loc) {
      loc[q] = v[k];  // no +1 here -- not ready to go to Matlab coords yet
    }
    out[q] = (q-v[k])*(q-v[k])+cost[v[k]];
  }
}

/****************************************************************************/
//
// helper function to do 2d minimum convolution
//
void gdt2d(double *cost, double *out, int *loc, int *v, double *z, 
            int nrow, int ncol, double x_multiply, double y_multiply, 
            double *scratch, double *scratch2, int *lscratch) {
  int i, j;
  double yx_ratio;
  int *loc2 = loc+(nrow*ncol);

  // pass over rows
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      scratch[j] = cost[i+j*nrow]*x_multiply;
    }
    if (loc) {
      gdt1d(scratch,scratch2+i*ncol,lscratch+i*ncol,v,z,ncol);
    } else {
      gdt1d(scratch,scratch2+i*ncol,0,v,z,ncol);
    }
  }

  // pass over columns
  yx_ratio = y_multiply/x_multiply;
  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
      scratch[i] = scratch2[i*ncol+j]*yx_ratio;
    }
    if (loc) {
      gdt1d(scratch,out+j*nrow,loc2+j*nrow,v,z,nrow);
    } else {
      gdt1d(scratch,out+j*nrow,0,v,z,nrow);
    }
  }
  for (i = 0; i < nrow*ncol; i++) {
    out[i] /= y_multiply;
  }

  // rewrite loc so that it contains actual indices.
  // original format has x referring to y.
  if (loc) {
    for (i = 0; i < nrow; i++) {
      for (j = 0; j < ncol; j++) {
        loc[i+j*nrow] = lscratch[((loc2[i+j*nrow]))*ncol+j];
      }
    }
  }
  //mexPrintf("Loc rewritten.\n");  
}


/****************************************************************************/
//
// helper function to do matrix translation --
// uses linear interpolation for subpixel accuracy
//
void tranDT(double *dt, double *dt2, int nrow, int ncol, 
            double x, double y, double maxd) {
  int i, j;
  int x1 = FLOOR(x);
  int y1 = FLOOR(y);
  int x2 = x1+1;
  int y2 = y1+1;
  int x3 = x1+2;
  int y3 = y1+2;
  int x0 = x1-1;
  int y0 = y1-1;
  double xf = x-x1;
  double yf = y-y1;
  double t0 = ((2-yf)*yf-1)*yf;
  double t1 = (3*yf-5)*yf*yf+2;
  double t2 = ((4-3*yf)*yf+1)*yf;
  double t3 = (yf-1)*yf*yf;
  double s0 = ((2-xf)*xf-1)*xf;
  double s1 = (3*xf-5)*xf*xf+2;
  double s2 = ((4-3*xf)*xf+1)*xf;
  double s3 = (xf-1)*xf*xf;
  double q00 = t0*s0/4;
  double q01 = t0*s1/4;
  double q02 = t0*s2/4;
  double q03 = t0*s3/4;
  double q10 = t1*s0/4;
  double q11 = t1*s1/4;
  double q12 = t1*s2/4;
  double q13 = t1*s3/4;
  double q20 = t2*s0/4;
  double q21 = t2*s1/4;
  double q22 = t2*s2/4;
  double q23 = t2*s3/4;
  double q30 = t3*s0/4;
  double q31 = t3*s1/4;
  double q32 = t3*s2/4;
  double q33 = t3*s3/4;

  for (i = 0; i < nrow; i++) {
    if ((i+y0 < 0)||(i+y3 >= nrow)) {
      for (j = 0; j < ncol; j++) {
        dt2[i+nrow*j] = maxd;
      }
    } else {
      for (j = 0; j < ncol; j++) {
        if ((j+x0 < 0)||(j+x3 >= ncol)) {
          dt2[i+nrow*j] = maxd;
        } else {
          dt2[i+nrow*j] = 
            q00*dt[i+y0+nrow*(j+x0)]+q01*dt[i+y0+nrow*(j+x1)]
            +q02*dt[i+y0+nrow*(j+x2)]+q03*dt[i+y0+nrow*(j+x3)]
            +q10*dt[i+y1+nrow*(j+x0)]+q11*dt[i+y1+nrow*(j+x1)]
            +q12*dt[i+y1+nrow*(j+x2)]+q13*dt[i+y1+nrow*(j+x3)]
            +q20*dt[i+y2+nrow*(j+x0)]+q21*dt[i+y2+nrow*(j+x1)]
            +q22*dt[i+y2+nrow*(j+x2)]+q23*dt[i+y2+nrow*(j+x3)]
            +q30*dt[i+y3+nrow*(j+x0)]+q31*dt[i+y3+nrow*(j+x1)]
            +q32*dt[i+y3+nrow*(j+x2)]+q33*dt[i+y3+nrow*(j+x3)];
        }
      }
    }
  }
}


/****************************************************************************/

// for debugging:
void printArr(double *arr, int m, int n) {
  for (int i = 0; i <m; i++) {
    for (int j = 0; j < n; j++) {
      mexPrintf("%g ",arr[i+j*m]);
    }
    mexPrintf("\n");
  }
}

/****************************************************************************/
//
// gateway driver to call the distance transform code from matlab
//
// This is the matlab entry point
void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i, j, k, c, nrow, ncol, maxdim, nmpt;
  int mxfn, myfn, mrfn, mvxfn, mvyfn, mpfn, mchfn;
  int root, curt, ch, prnt, dtarrbase, dtarrlim;
  int ldim[3], *v, *mp, *mnch, *mdone, *mndesc, *dtarrhist;
  int *loc, *loc2, *lscratch, **locarr, **locarr2;
  char *pbimgc;
  double x_multiply, y_multiply, yx_ratio, maxd, *z;
  double *cost, *out, *scratch, *scratch2, *scl;
  double *mx, *my, *mr, *mvx, *mvy, *pbimg, *dt1, *tranarr;
  double **mch, **locp, **dtarr, **dtarrp;
  mxArray *cell, *field;

  // check for proper number and size of arguments
  errCheck(nrhs >= 2,"Arguments:  model, bimg, [maxd]");
  errCheck(nrhs <= 3,"Arguments:  model, bimg, [maxd]");
  errCheck(nlhs <= 3,"Outputs:  dt, [loc], [scale]");
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
  if (nrhs>2) {
    errCheck(mxIsDouble(prhs[2])&&!mxIsComplex(prhs[2]),
             "maxd must be real double.");
    errCheck(mxGetNumberOfElements(prhs[2])==1,
             "maxd must be scalar.");
    maxd = *mxGetPr(prhs[2]);
  } else {
    maxd = LARGE;
  }
  nmpt = (int)mxGetNumberOfElements(prhs[0]);
  nrow = (int)mxGetM(prhs[1]);
  ncol = (int)mxGetN(prhs[1]);

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
    loc2 = 0;
    locp = 0;
  }

  // transcribe image
  pbimg = mxGetPr(prhs[1]);
  pbimgc = (char*)pbimg;
  dt1 = out;  // this is where the answer will end up, eventually
  if (mxIsDouble(prhs[1])) {
    for (int i = 0; i <nrow*ncol; i++) {
      if (pbimg[i]) {
        dt1[i] = 0;
      } else {
        dt1[i] = maxd;
      }
    }
  } else {
    for (int i = 0; i <nrow*ncol; i++) {
      if (pbimgc[i]) {
        dt1[i] = 0;
      } else {
        dt1[i] = maxd;
      }
    }
  }

  // allocate scratch space for distance transform
  maxdim = MAX(nrow,ncol);
  scratch = (double*)mxMalloc(maxdim*sizeof(double));
  scratch2 = (double*)mxMalloc(nrow*ncol*sizeof(double));
  if (loc) {
    lscratch = (int*)mxMalloc(nrow*ncol*sizeof(int));
  } else {
    lscratch = 0;
  }
  v = (int*)mxMalloc(maxdim*sizeof(int));
  z = (double*)mxMalloc((maxdim+1)*sizeof(double));

  // allocate other variables
  mdone = (int*)mxMalloc(nmpt*sizeof(int));
  mndesc = (int*)mxMalloc(nmpt*sizeof(int));
  tranarr = (double*)mxMalloc(nrow*ncol*sizeof(double));  // stores translation
  dtarr = (double**)mxMalloc(nmpt*sizeof(double*));  // stores intermediate dt
  dtarrp = (double**)mxMalloc(nmpt*sizeof(double*));  // pointer to storage
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
      locarr[i] = (int*)mxMalloc(2*nrow*ncol*sizeof(int));
      locarr2[i] = locarr[i]+nrow*ncol;
    }
  } else {
    for (i = 0; i < nmpt; i++) {
      locarr[i] = 0;
      locarr2[i] = 0;
    }
  }

  // compute distance transform on input image
  gdt2d(dt1,dt1,0,v,z,nrow,ncol,1,1,scratch,scratch2,lscratch);

  // run computation
  curt = root;
  dtarrlim = 0;
  dtarrbase = 0;
  dtarrp[root] = dt1;
  k = 0;
  while (curt != -1) {

    // set base point
    if (mdone[curt]==0) {
      dtarrhist[curt] = dtarrbase;
    }

    // figure out if we need to descend to children
    if (mdone[curt]<mnch[curt]) {
      // figure out current child
      ch = (int)(mch[curt][mdone[curt]])-1;  // Matlab indices start from 1
     
      // we haven't yet computed the current child; set up to do that now
      ch = (int)(mch[curt][mdone[curt]])-1;  // Matlab indices start from 1
      
      // allocate new dt storage space if necessary
      if (dtarrlim==dtarrbase) {
        dtarr[dtarrlim] = (double*)mxMalloc(nrow*ncol*sizeof(double));
        dtarrlim++;
      }
      
      // copy dt1
      for (i = 0; i < nrow*ncol; i++) {
        dtarr[dtarrbase][i] = dt1[i];
      }
           
      // set up storage for this child
      dtarrp[ch] = dtarr[dtarrbase];
      dtarrbase++;
      
      // move on, but when we return, remember that this child is done
      mdone[curt]++;
      curt = ch;
    } else {
      // all children complete; assemble result
      for (c = 0; c < mnch[curt]; c++) {
        ch = (int)(mch[curt][c])-1;  // Matlab indices start from 1
        errCheck(dtarrp[ch]!=0,"Null pointer");
        mxAssert((ch>=0)&&(ch<nmpt),"Out of bounds");
        tranDT(dtarrp[ch],tranarr,nrow,ncol,mx[ch],my[ch],maxd*mndesc[ch]);
        gdt2d(tranarr,dtarrp[ch],locarr[ch],v,z,nrow,ncol,mvx[ch],mvy[ch],
               scratch,scratch2,lscratch);
        for (i = 0; i < nrow*ncol; i++) {
          dtarrp[curt][i] += dtarrp[ch][i];
        }
      }

      // now go back up to the parent
      if (mp[curt]>0) {
        mxAssert((curt>=0)&&(curt<nmpt),"Out of bounds");
        mxAssert((mp[curt]-1>=0)&&(mp[curt]-1<nmpt),"Out of bounds");
        mndesc[mp[curt]-1] += mndesc[curt];
      }
      dtarrbase = dtarrhist[curt];
      curt = mp[curt]-1;  // Matlab indices start from 1
    }
  }

  // do final translation and distance transform
  tranDT(dt1,tranarr,nrow,ncol,mx[root],my[root],maxd*nmpt);
  gdt2d(tranarr,dt1,locarr[root],v,z,nrow,ncol,mvx[root],mvy[root],
         scratch,scratch2,lscratch);

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
        // copy results for this node
        for (i = 0; i < nrow; i++) {
          for (j = 0; j < ncol; j++) {
            if (prnt == -1) {
              locp[i+j*nrow][2*curt] = BOUND(j+x,0,ncol-1);
              locp[i+j*nrow][2*curt+1] = BOUND(i+y,0,nrow-1);
            } else {
              int jj = (int)BOUND(locp[i+j*nrow][2*prnt],0,ncol-1);
              int ii = (int)BOUND(locp[i+j*nrow][2*prnt+1],0,nrow-1);
              locp[i+j*nrow][2*curt] = locarr[curt][ii+nrow*jj]+x;
              locp[i+j*nrow][2*curt+1] = locarr2[curt][ii+nrow*jj]+y;
            }
          }
        }
      }

      // descend to children if any left
      if (mdone[curt]<mnch[curt]) {
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
  }

  // compute scale if needed
  if (nlhs>=3) {
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
    mxFree(dtarr[i]);
  }
  if (loc) {
    for (i = 0; i < nmpt; i++) {
      mxFree(locarr[i]);
    }
  }
  mxFree(locarr);
  mxFree(dtarrhist);
  mxFree(dtarrp);
  mxFree(dtarr);
  mxFree(mndesc);
  mxFree(mdone);
  mxFree(z);
  mxFree(v);
  if (locp) {
    mxFree(lscratch);
  }
  mxFree(scratch2);
  mxFree(scratch);
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
}
