#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WIN 1000 // window size

/*******************************************************************/
int state = 0;      // tracer state idle/working - must be 0 in the beginning
double f0, a0, p0;  // initial conditions
double t0 = 0;      // for absolute time output
int slope;          // signal slope (+1/-1)
double xf1, xf2;    // filtering oscillator state

FILE *flt_dbg_file = NULL; // for filtered signal and filter state
FILE *ph_dbg_file  = NULL; // calculated signal before adjusting phase and phase diff
FILE *amp_dbg_file = NULL; // calculated signal before adjusting amp and amp diff
FILE *fin_dbg_file = NULL; // final calculated signal

struct result_t{
  double time, amp, freq, err;
} res;


/*******************************************************************/
void
process_window(double *xf, int np, double dt){
  double level=0.5;   // some amplitude level
  int j;
  double aS=0, aSx=0, aSy=0, aSxx=0, aSxy=0;
  double pS=0, pSx=0, pSxx=0, pSxxx=0, pSxxxx=0, pSy=0, pSxy=0, pSxxy=0;

  double pdet, dp0,dp1,dp2,da0,da1;
  double err; // step error

  double slope_lock; // parameter for locking slope at short times
  double period;     // signal period


  /*******************************************************************/
  // idle -> run, calculate initial values
  if (state==0){
    double x0, x1;
    double base;
    double bS=0, bSy=0;
    int iper;

    // find mean value
    for (j=0; j<np; j++) {bSy+=xf[j]; bS++; }
    base=bSy/bS;

    // find amplitude
    a0 = 0;
    for (j=0; j<np; j++){
      if (fabs(xf[j]-base)>a0) a0 = fabs(xf[j]-base);
    }

    x0=(xf[0]-base)/a0;
    x1=(xf[1]-base)/a0;
    // find initial frequency
    period=0;
    slope=( (x1>x0 || x0<-level) && x0<level)? 1:-1; // initial slope
    for (j=0; j<np; j++){
      double x=(xf[j]-base)/a0;
      if (slope>0 && x>level)  { slope = -1; }
      if (slope<0 && x<-level) { slope = 1; period+=1; }
    }
    if (period<1) period=1;
    f0 = period/np/dt;
    p0 = 0;

    // set initial conditions for filtering oscillator (phase shift!)
    iper = 1/f0/dt; // points/period
    if (iper<4 || iper>np) {t0 += j*dt; return;}
    xf1=xf[3*iper/4-1]/a0;
    xf2=xf[3*iper/4]/a0;

    // go to state 1 (run!)
  }
  state++;

  /*******************************************************************/
  // first loop - filtering + phase adjastment
  if (flt_dbg_file)
    fprintf(flt_dbg_file, "\n#%13s %14s %14s\n",
      "t", "x", "dx/dt", "damping");
  if (ph_dbg_file)
    fprintf(ph_dbg_file, "\n#%13s %14s %14s %14s\n",
      "t", "approx sig", "phase diff", "slope+level");
  period=0;
  slope_lock=0;


  for (j=0; j<np; j++){
    double t = j*dt;
    double w = 2*M_PI*f0;
    double pp = p0 + w*t;
    double xp = a0*sin(pp); // approximate signal
    double xma, damp;

    // Filtering oscillator with frequency f0 and damping damp.
    // In the beginning of the signal damping is 1 (for filter
    // stabilization) and after ~1 period it drops
    // to 2/(number of periods per window)
    damp=0.5/f0/np/dt;
    if (state==1) damp += exp(-1.0*j*f0*dt);

    double xf3 = (xf2*(2-dt*dt*w*w) - xf1*(1-dt*w*damp) +
            2*damp* xf[j]/a0*dt*dt*w*w)/(1+dt*w*damp);
    // note: driving firce must be smaller at low damping to
    // keep amplitude constant.
    xf[j] = a0*(xf3-xf1)/2/dt/w; // oscillator velocity
    if (flt_dbg_file)
      fprintf(flt_dbg_file, "%14e %14e %14e %14e\n",
        t0+j*dt, xf2*a0, xf[j], damp);
    xf1=xf2; xf2=xf3;

    xma = xf[j]/a0; // 'real' sin(ph)

    // switch slope if needed
    if (slope>0 && (
         (xma>level && t-slope_lock>1/(8*f0)) ||
         (t-slope_lock>1/(1.5*f0)) )) {
      slope = -1; slope_lock=t;
    }
    if (slope<0 && (
         (xma<-level && t-slope_lock>1/(8*f0)) ||
         (t-slope_lock>1/(1.5*f0)) )) {
      slope = 1; slope_lock=t; period+=1;
    }

    // calculate phase difference
    if (fabs(xma)<level) {
      double y = 2*M_PI*period + (slope>0 ? asin(xma):(M_PI-asin(xma))) - pp; // dp
      if (!isnan(y)){ pS++; pSx+=t; pSxx+=t*t; pSxxx+=t*t*t; pSxxxx+=t*t*t*t;
                      pSy+=y; pSxy+=t*y; pSxxy+=t*t*y; }
      if (ph_dbg_file)
        fprintf(ph_dbg_file, "%14e %14e %14e %14e\n",
          t0+j*dt, xp, y, slope*a0*level);
    }
  }

  // Phase adjastments:
  // quadratic fit: dph(t) = dph0 + dph1*t + dph2*t*t
  pdet = (pSxxxx*pSxx*pS + 2*pSxxx*pSxx*pSx
         - pSxx*pSxx*pSxx - pSxxx*pSxxx*pS - pSx*pSx*pSxxxx);
  dp0 =  (pSxxxx*pSxx*pSy + pSxxx*pSxx*pSxy + pSxxx*pSxxy*pSx
         - pSxx*pSxx*pSxxy - pSxxx*pSxxx*pSy - pSxy*pSx*pSxxxx)/pdet;
  dp1 = -(pSxxxx*pSx*pSy + pSxx*pSxx*pSxy + pSxxx*pSxxy*pS
         - pSxxy*pSxx*pSx - pSxxx*pSxx*pSy - pSxy*pS*pSxxxx)/pdet;
  dp2 =  (pSxxx*pSx*pSy + pSxx*pSxy*pSx + pSxxy*pSxx*pS
         - pSxxy*pSx*pSx - pSxx*pSxx*pSy - pS*pSxy*pSxxx)/pdet;


  /*******************************************************************/
  // second loop - amplitude adjastment
  if (amp_dbg_file)
    fprintf(amp_dbg_file, "\n#%13s %14s %14s\n",
      "t", "approx sig", "amp diff");
  for (j=0; j<np; j++){
    double t = j*dt;
    double pp = p0 + 2*M_PI*f0*t + dp0 + dp1*t + dp2*t*t;
    double xp = a0*sin(pp);

    if (fabs(xp/a0) >= level) {
      double y = (xf[j]-xp)/sin(pp); // da
      aS++; aSx+=t; aSy+=y; aSxx+=t*t; aSxy+=t*y;
      if (amp_dbg_file)
        fprintf(amp_dbg_file, "%14e %14e %14e\n",
          t0+j*dt, xp, y);
    }
  }

  // Amplitude adjastments:
  if (state==1){
    // linear fit da(t) = da0 + da1(t) -- good for first step
    da0 = (aSxy*aSx-aSxx*aSy) / (aSx*aSx-aSxx*aS);
    da1 = (aSx*aSy-aSxy*aS) / (aSx*aSx-aSxx*aS);
  }
  else {
    // linear fit da(t) = da1(t) -- if initial amp is good
    da0 = 0;
    da1 = aSxy/aSxx;
  }
  //  // da(t) = da0 -- also possible
  //  da0 = aSy/aS;
  //  da1 = 0;

  /*******************************************************************/
  // calculate error and print filtered signal
  err=0; aS=0;
  if (fin_dbg_file)
    fprintf(fin_dbg_file, "\n#%13s %14s\n", "t", "approx sig");
  for (j=0; j<np; j++){
    double t = j*dt;
    double pp = p0 + 2*M_PI*t*f0 + dp0 + dp1*t + dp2*t*t;
    double xp = (a0+da0+da1*t)*sin(pp);
    err+=(xf[j]-xp) * (xf[j]-xp); aS++;
    if (fin_dbg_file)
      fprintf(fin_dbg_file, "%14e %14e\n", t0+j*dt, xp);
  }

  err = sqrt(err/aS);
  if (err > a0+da0 + 0.5*j*dt * da1 ||
      isnan(da0) || isnan(da1) ||
      isnan(dp0) || isnan(dp1) || isnan(dp2)) {state=0; t0 += j*dt;
//fprintf(stderr,">> %e %e %e %e %e\n", da0, da1, dp0, dp1, dp2);
       return;
  }

  // output data (middle of the window)
  res.time  = t0 + 0.5*j*dt;
  res.freq  = f0 + (dp1 + dp2*j*dt)/2/M_PI;
  res.amp   = a0 + da0 + 0.5*j*dt * da1;
  res.err   = err;

  // initial values for the next step
  p0 += dp0 + j*dt * dp1  + dp2*j*dt*j*dt;
  f0 += (dp1 + 2*dp2*j*dt)/2/M_PI;
  a0 += da0 + j*dt * da1;
  t0 += j*dt;
}




/*******************************************************************/
int
main(){
  int i, fmt, val;
  double dt, dx;
  double xm[WIN]; // initial signal


//  flt_dbg_file = fopen("data_flt.dat", "w");
//  ph_dbg_file  = fopen("data_ph.dat", "w");
//  amp_dbg_file = fopen("data_amp.dat", "w");
  fin_dbg_file = fopen("data_fin.dat", "w");

  FILE * tmp_file = fopen("data_tmp.dat", "w");


  { // read signal header
    char *line, *str, *tok, *saveptr;
    int j, len=0;

    if (getline(&line, &len, stdin)<1) exit(1);
    fprintf(stderr, ">>> %s\n", line);
    for (j=1,str=line; ;j++,str=NULL){
      tok = strtok_r(str, ",", &saveptr);
      if (tok == NULL) break;
      switch (j){
        case 1: fmt=atoi(tok); break;
        case 5: dt=atof(tok); break;
        case 6: t0=atof(tok); break;
        case 8: dx=atof(tok); break;
      }
    }
    // skip two lines
    for (j=0;j<2;j++) getline(&line, &len, stdin);
    free(line);
  }

  printf("#%13s %14s %14s %14s\n",
    "time", "freq", "amp", "err");
  do {
    for (i=0; i<WIN; i++) {
      scanf("%i", &val);
      if (fmt==0) val=(val>=0)?val-127:val+127;
      xm[i] = val*dx;
      fprintf(tmp_file, "%d\n", val);

      if (feof(stdin)) break;
    }
    process_window(xm, i, dt);
    fprintf(stderr, state>0?"o":".");

    if (state)
      printf("%14e %14e %14e %14e\n",
        res.time, res.freq, res.amp, res.err);

  } while (!feof(stdin));

  return 0;
}
