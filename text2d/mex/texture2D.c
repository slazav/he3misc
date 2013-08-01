#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  int m_in, n_in, size;
  int npt=0, nptr=0, ns=0;
  double *textpar=NULL, *specpar=NULL;
  double *textur=NULL, *resspec=NULL;
  double *apsipar=NULL;
  int initype=1, msglev=-10;
  void calctexture_(int*,double*,int*,double*,int*,double*,double*,int*,double*,int*);
  
  if (nrhs < 6)
    mexErrMsgTxt("7 inputs required.");
  /* 1st input argument is the number of discretization points
     for texture, or empty (then do not optimize the texture,
     use initial value)

     2nd input argument is texture calculation parameters
     1 - temperature / Tc
     2 - pressure, bar
     3 - larmor frequency, kHz
     4 - cylinder radius, cm
     5 - rotation velocity, rad/s
     6 - omega_v, rad/s
     7 - lambda/omega (s/rad) if >=0
         use calculated lambda/omega if == -1
     8 - lambga_HV (kg/(m^3 T^2)) if >= 0
         use calculated lambga_HV if == -1
     9 - ksi (dimensionless,same scale as in the func that is used to calculate ksi)

     3rd input argument is the number of points in the output spectrum,
     or empty if spectrum is not needed

     4th input argument is spectrum calculation parameters
     1 - half-width of NMR line
     2 - margin for automatic region determination

     5th input argument is initial condition
     1 - texture w/o 90 deg peak
     2 - texture with 90 deg peak
     n by 3 array - 1st column ignored , 2nd - alpha, 3rd - beta

     6th input argument is the sqrt(2)*sin(B_µ/2)=(A*|Psi|) -vector

     1st output argument
     texture - n by 3 array - 1st column r, 2nd - alpha, 3rd - beta, 4th energies
     
     2nd output argument
     spectrum - m by 2 array - 1st column freq, 2nd column absorption
  */
  
  m_in = mxGetM(prhs[0]);
  n_in = mxGetN(prhs[0]);
  size = m_in * n_in;
  if (size == 0) npt = 0;
  else {
    if (size != 1 || !mxIsDouble(prhs[0]))
      mexErrMsgTxt("1st input must be a scalar.");
    else npt = *(double*)mxGetPr(prhs[0]);
  }
  
  m_in = mxGetM(prhs[1]);
  n_in = mxGetN(prhs[1]);
  size = m_in * n_in;
  if (size != 10 || !mxIsDouble(prhs[1]))
    mexErrMsgTxt("2nd input must be double vector of length 10.");
  else textpar = (double*)mxGetPr(prhs[1]);
  
  m_in = mxGetM(prhs[2]);
  n_in = mxGetN(prhs[2]);
  size = m_in * n_in;
  if (size == 0) ns = 0;
  else {
    if (size != 1 || !mxIsDouble(prhs[2]))
      mexErrMsgTxt("3rd input must be a scalar.");
    else ns = *(double*)mxGetPr(prhs[2]);
  }
  
  if (ns > 0) {
    m_in = mxGetM(prhs[3]);
    n_in = mxGetN(prhs[3]);
    size = m_in * n_in;
    if (size != 2 || !mxIsDouble(prhs[3]))
      mexErrMsgTxt("4th input must be double vector of length 2.");
    else specpar = (double*)mxGetPr(prhs[3]);
  }

  if (!mxIsDouble(prhs[4]))
      mexErrMsgTxt("5th input must be double array.");
  m_in = mxGetM(prhs[4]);
  n_in = mxGetN(prhs[4]);
  if (m_in == 1 && n_in == 1 && npt > 0)
    initype = *(double*)mxGetPr(prhs[4]);
  else {
    if (n_in != 3 || (npt && npt != m_in -1))
      mexErrMsgTxt("Wrong initial condition.");
    if (npt) initype = 3;
    else {
      npt = m_in - 1;
      initype = 4;
    }
  }

  plhs[0] = mxCreateDoubleMatrix(npt+1, 4, mxREAL);
  textur = mxGetPr(plhs[0]);
  if (initype >= 3)
    memcpy(textur,mxGetPr(prhs[4]),(npt+1)*3*sizeof(double));
  if (nlhs > 1 && ns > 0) {
    plhs[1] = mxCreateDoubleMatrix(ns+1, 2, mxREAL);
    resspec = mxGetPr(plhs[1]);
  } else
    ns = 0;


  m_in = mxGetM(prhs[5]);
  n_in = mxGetN(prhs[5]);
  size = m_in * n_in;
  if (size !=npt+1  || !mxIsDouble(prhs[5]))
    mexErrMsgTxt("6th input must be double vector of length N+1.");
  else apsipar = (double*)mxGetPr(prhs[5]);


  m_in = mxGetM(prhs[6]);
  n_in = mxGetN(prhs[6]);
  size = m_in * n_in;
  if (size == 0) nptr = 0;
  else {
    if (size != 1 || !mxIsDouble(prhs[6]))
      mexErrMsgTxt("7th input must be a scalar.");
    else nptr = *(double*)mxGetPr(prhs[6]);
  }

  calctexture_(&npt,textpar,&ns,specpar,&initype,textur,resspec,&msglev,apsipar,&nptr);

  if (*textur == -1)
    mexErrMsgTxt("Too big number of discretization intervals.");
  else if (*textur == -2)
    mexErrMsgTxt("Error in interpolation, aborting.");
}

