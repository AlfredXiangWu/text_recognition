/****************************************************************************/
//
// Matlab interface file:  gdt2d.cpp
//
// Written 6/10 by N. Howe, inspired by Matlab code by R. Fergus.
// Computes 2-d generalized distance transform.
//
/****************************************************************************/

#include "mex.h"
#include "mymex.h"

#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree

// macro for computing maximum
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

// macro for checking error conditions
#ifndef errCheck
#define errCheck(a,b) if (!(a)) mexErrMsgTxt((b));
#endif

/****************************************************************************/
//
// helper function to do 1d minimum convolution
//
void psdt1d(double *cost, double *out, double *loc, int *v, double *z, int size_x) {
  int k, q;
  double s;

  // set up
  k = 0;
  v[0] = 0;
  z[0] = -mxGetInf();
  z[1] = mxGetInf();

  // compute
  for (q=1; q<size_x; q++) {
    mxAssert(q<size_x,"Q out of range");
    mxAssert(k<size_x,"v out of range");
    mxAssert(v[k]+1<size_x,"cost out of range");
    s = ((cost[q]+q*q)-(cost[v[k]]+v[k]*v[k]))/(2*(q-v[k]));  // intercept
    while (s<=z[k]) {
      k = k-1;
      mxAssert(q<size_x,"Q out of range");
      mxAssert(k<size_x,"v out of range");
      mxAssert(v[k]+1<size_x,"cost out of range");
      s = ((cost[q]+q*q)-(cost[v[k]]+v[k]*v[k]))/(2*(q-v[k]));
    }
    k = k+1;
    mxAssert(k<size_x,"k out of range");
    v[k] = q;
    z[k] = s;
    z[k+1] = mxGetInf();
  }
  k = 0;
  for (q=0; q<size_x; q++) {
    while (z[k+1]<q) {
      k++;
    }
    if (loc) {
      loc[q] = v[k]+1;
    }
    out[q] = (q-v[k])*(q-v[k])+cost[v[k]];
  }
}

/****************************************************************************/
//
// gateway driver to call the distance transform code from matlab
//
// This is the matlab entry point
void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i, j, nrow, ncol, maxdim;
  int ldim[3], *v;
  double x_multiply, y_multiply, yx_ratio, *z;
  double *cost, *out, *loc, *loc2, *scratch, *lscratch, *scratch2;

  // check for proper number and size of arguments
  errCheck(nrhs >= 1,"Arguments:  cost_function, [x_multiply], [y_multiply]");
  errCheck(nrhs <= 3,"Arguments:  cost_function, [x_multiply], [y_multiply]");
  errCheck(mxIsDouble(prhs[0]),"cost_function must be double.");
  errCheck(!mxIsComplex(prhs[0]),"cost_function must be real.");
  if (nrhs>=2) {
    errCheck(mxIsDouble(prhs[1]),"x_multiply must be double.");
    errCheck(!mxIsComplex(prhs[1]),"x_multiply must be real.");
    errCheck(mxGetNumberOfElements(prhs[1])==1,"x_multiply must be scalar");
    x_multiply = *mxGetPr(prhs[1]);
  } else {
    x_multiply = 1;
  }
  if (nrhs>=3) {
    errCheck(mxIsDouble(prhs[2]),"y_multiply must be double.");
    errCheck(!mxIsComplex(prhs[2]),"y_multiply must be real.");
    errCheck(mxGetNumberOfElements(prhs[2])==1,"y_multiply must be scalar");
    y_multiply = *mxGetPr(prhs[2]);
  } else {
    y_multiply = 1;
  }
  errCheck(nlhs <= 2,"Outputs:  out, [loc]");
  errCheck(mxGetNumberOfDimensions(prhs[0]) <= 2,"cost_function must be 2D");
  nrow = (int)mxGetM(prhs[0]);
  ncol = (int)mxGetN(prhs[0]);

  // set up input variables
  cost = mxGetPr(prhs[0]);

  // allocate output space
  plhs[0] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  out = mxGetPr(plhs[0]);
  if (nlhs==2) {
    ldim[0] = nrow;
    ldim[1] = ncol;
    ldim[2] = 2;
    plhs[1] = mxCreateNumericArray(3, ldim, mxDOUBLE_CLASS, mxREAL);
    loc = mxGetPr(plhs[1]);
    loc2 = loc+(nrow*ncol);
  } else {
    loc = loc2 = 0;
  }

  // allocate scratch space
  maxdim = MAX(nrow,ncol);
  scratch = (double*)mxMalloc(maxdim*sizeof(double));
  scratch2 = (double*)mxMalloc(nrow*ncol*sizeof(double));
  if (loc) {
    lscratch = (double*)mxMalloc(nrow*ncol*sizeof(double));
  } else {
    lscratch = 0;
  }
  v = (int*)mxMalloc(maxdim*sizeof(int));
  z = (double*)mxMalloc((maxdim+1)*sizeof(double));

  // pass over rows
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      scratch[j] = cost[i+j*nrow]*x_multiply;
    }
    if (loc) {
      psdt1d(scratch,scratch2+i*ncol,lscratch+i*ncol,v,z,ncol);
    } else {
      psdt1d(scratch,scratch2+i*ncol,0,v,z,ncol);
    }
  }

  // pass over columns
  yx_ratio = y_multiply/x_multiply;
  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
      scratch[i] = scratch2[i*ncol+j]*yx_ratio;
    }
    if (loc) {
      psdt1d(scratch,out+j*nrow,loc2+j*nrow,v,z,nrow);
    } else {
      psdt1d(scratch,out+j*nrow,0,v,z,nrow);
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
        loc[i+j*nrow] = lscratch[(int)((loc2[i+j*nrow])-1)*ncol+j];
      }
    }
  }

  // free stuff
  mxFree(z);
  mxFree(v);
  if (loc) {
    mxFree(lscratch);
  }
  mxFree(scratch2);
  mxFree(scratch);
}
