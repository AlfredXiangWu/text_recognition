/****************************************************************************/
/*
 * Matlab C routine header file:  mymex
 *
 * Written 6/99 by N. Howe.
 */
/****************************************************************************/

#ifndef MYMEX_H
#define MYMEX_H

#include <mex.h>

#define malloc mxMalloc
#define calloc mxCalloc
#define free mxFree

void store_field(mxArray *target, int id, int fnum, double val);
#define mxIsScalar(x) (mxIsNumeric(x)&&(mxGetM(x)==1)&&(mxGetN(x)==1))

#define errCheck(a,b) if (!(a)) mexErrMsgTxt((b));

#define mxGetFieldByName(m,i,f) (mxGetFieldByNumber((m),(i),mxGetFieldNumber((m),(f))))

#endif
