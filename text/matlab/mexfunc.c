#include "mex.h"
#include <string.h>
#include <stdio.h>

#define MAXN 10
struct text_pars{
  double alpha[MAXN], beta[MAXN];
  double a,d;
};

/* Empty constructor for texture parameter structure.
   Set arrays to 0 and all other values to -1.
   Move to the library! */
void
text_init_zero(struct text_pars *cst){
  memset(cst->alpha, 0, sizeof(cst->alpha));
  memset(cst->beta,  0, sizeof(cst->beta));
  cst->a=-1.0;
  cst->d=-1.0;
}

/**********************************************************/
/* set a field in matlab/octave structure from an array */
void
set_field(mxArray *mst, int nfld, int N, double *arr){
  int i;
  double *arr1 = mxMalloc(N*sizeof(double));
  mxArray *mxarr = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
  for (i=0; i<N; i++) arr1[i]=arr[i];
  mxSetData(mxarr, arr1);
  mxSetFieldByNumber(mst, 0, nfld, mxarr);
}

/* get a field in matlab structure into array */
void
get_field(const mxArray *mst, int nfld, int N, double *arr){
  int i;
  mxArray *v;
  double *arr1;
  v = mxGetFieldByNumber(mst, 0, nfld);
  if (!mxIsDouble(v))
    mexErrMsgTxt("Not a floating point array in a structure field");
  if (mxGetNumberOfElements(v)!=N)
    mexErrMsgTxt("Wrong size of array in a structure field");
  arr1 = mxGetData(v);
  for (i=0;i<N;i++) arr[i] = arr1[i];
}

/**********************************************************/
/* convert c structure to matlab */
mxArray *
pars_c2mat(struct text_pars *cst){
  int nkeys = 4;
  const char *keys[] = {"alpha", "beta", "a", "d"};
  mxArray * mst = mxCreateStructMatrix(1,1,nkeys, keys);

  set_field(mst, 0, MAXN, cst->alpha);
  set_field(mst, 1, MAXN, cst->beta);
  set_field(mst, 2, 1, &cst->a);
  set_field(mst, 3, 1, &cst->d);
  return mst;
}

/* convert matlab structure to c */
struct text_pars
pars_mat2c(const mxArray *mst){
  int i;
  const char * k;
  struct text_pars cst;

  if (!mxIsStruct(mst))
    mexErrMsgTxt("Structure expected");
  if (mxGetNumberOfElements(mst)!=1)
    mexErrMsgTxt("One element in StructureArray is needed");

  text_init_zero(&cst);

  for (i = 0; i < mxGetNumberOfFields(mst); i++){
    k = mxGetFieldNameByNumber(mst, i);
    if (strcmp(k, "alpha")==0){ get_field(mst, i, MAXN, cst.alpha); continue;}
    if (strcmp(k, "beta")==0) { get_field(mst, i, MAXN, cst.beta);  continue;}
    if (strcmp(k, "a")==0)    { get_field(mst, i, 1, &cst.a);  continue;}
    if (strcmp(k, "d")==0)    { get_field(mst, i, 1, &cst.d);  continue;}
    mexErrMsgTxt("Unknown field");
  }
  return cst;
}

/**********************************************************/
void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){

  struct text_pars pars;

  if ((nrhs != 1) || (!mxIsStruct(prhs[0])))
    mexErrMsgTxt("One input argument is needed");
  if (nlhs != 1)
    mexErrMsgTxt("One output argument is needed");

  pars = pars_mat2c(prhs[0]);

  plhs[0] = pars_c2mat(&pars);

  return;
}
