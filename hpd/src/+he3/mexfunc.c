#include "mex.h"
#include "math.h"
#include "../he3.h"


void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){
  int m_in, n_in, size, i;
  int npt=0, ns=0;
  double *in, *out;

/* Constants */
#if NARGIN == 0
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  out = (double*)mxGetPr(plhs[0]);
  *out = he3_const_.FUNC;
  return;
#endif

/* Functions of one argument */
#if NARGIN == 1
  if (nrhs != 1)
    mexErrMsgTxt("1 input required."); 

  m_in = mxGetM(prhs[0]);
  n_in = mxGetN(prhs[0]);
  size = m_in * n_in;
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("input must be a vector of double.");

  in = (double*)mxGetPr(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m_in, n_in, mxREAL);
  out = (double*)mxGetPr(plhs[0]);

  for (i=0; i<size; i++){
    out[i] = FUNC(in+i);
    if (out[i]<0) out[i]=NAN;
  }
  return;
#endif
}
