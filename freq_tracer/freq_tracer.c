#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// number of points
#define PTS 50000
#define WIN 1000
int
main(){
  int i;
  double dt=1.5e-6;   // time step
  int state = 0;      // tracer state idle/working - must be 0 in the beginning
  double np = 1000;   // window (number of points)
  double level=0.5;   // some amplitude level

  double xm[PTS]; // initial signal
  double xf[WIN]; // filtered signal - one window

  // initial conditions
  double f0, a0, p0;
  // step error
  double err;


  int filter=1; // do freq filtering
  double xf1, xf2; // filtering oscillator state

  int print_data = 1; // print debugging information to data2.dat
                      // 0  - nothing
                      // 1  - filtered signal
                      // 2  - calculated signal before adjusting phase
                      // 21 - phase adjastment
                      // 22 - period tracking during phase adjastment
                      // 3  - calculated signal before adjusting amp
                      // 4  - final calculated signal

  // phase locking parameters
  int slope;         // signal slope
  double slope_lock; // parameter for locking slope at short times
  double period;     // signal period

  FILE *F1 = fopen("data1.dat", "w");
  FILE *F2 = fopen("data2.dat", "w");
  FILE *F3 = fopen("data3.dat", "w");

  if (!F1 || !F2 || !F3) return 1;

  // create model signal
  {
    double amp0=1.2;     // signal amplitude
    double amp0n = 0.22; // noise amplitude
    double fre0=34567.0; // final frequency
    double df0=1234.0;   // frequency change
    double ta=0.012;     // amplitude relaxation time
    double tf=0.021;     // frequency relaxation time
    double ph=0.4;       // initial phase
    double t0=0.0011;    // start of signal
    double base=0;

    for (i=0; i<PTS; i++){
      double t=i*dt;
      double noise = amp0n * (double)random() / RAND_MAX - amp0n/2;
      double fre = fre0 + df0 * exp(-t/tf);
      double amp = t<t0? 0 : amp0 * exp(-t/ta);
      if (i>0) ph=ph+2*M_PI*fre*dt;
      while (ph>M_PI) ph-=2*M_PI;
      xm[i] = noise + amp * sin(ph) + base;
      fprintf(F1, "%14e %14e %14e %14e %14e\n", t, xm[i], fre, amp, amp0n);
    }
  }

  for (i=0; i<PTS; i++){
    int j;
    int np=WIN;

    double aS=0, aSx=0, aSy=0, aSxx=0, aSxy=0;
    double pS=0, pSx=0, pSxx=0, pSxxx=0, pSxxxx=0, pSy=0, pSxy=0, pSxxy=0;

    double pdet, dp0,dp1,dp2,da0,da1;
    double base;

    // copy part of the signal to xf
    if (i+np > PTS) np=PTS-i;
    for (j=i; j<i+np; j++) xf[j-i]=xm[j];

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
        if (j>=PTS) break;
        if (slope>0 && x>level)  { slope = -1; }
        if (slope<0 && x<-level) { slope = 1; period+=1; }
      }
      if (period<1) period=1;
      f0 = period/np/dt;
      p0 = 0;

      // set initial conditions for filtering oscillator (phase shift!)
      iper = 1/f0/dt; // points/period
      if (iper<4) iper=4;
      if (iper>np) iper=np;
      xf1=xf[3*iper/4-1]/a0;
      xf2=xf[3*iper/4]/a0;

      // go to state 1 (run!)
      state = 1;
    }

    period=0;
    slope_lock=0;
    // first loop - adjast phase
    for (j=0; j<np; j++){
      double t = j*dt;
      double w = 2*M_PI*f0;
      double pp = p0 + w*t;
      double xp = a0*sin(pp);
      double xma;

      // filtering oscillator
      if (filter){
        double g=0.5;
        double xf3 = (xf2*(2-dt*dt*w*w) - xf1*(1-dt*w*g) +
                 xf[j]/a0*dt*dt*w*w)/(1+dt*w*g);
        xf[j] = a0*(xf3-xf1)/2/dt/w; // oscillator velocity!
        if (print_data==1) fprintf(F2, "%14e %14e %14e\n", (i+j)*dt, xf[j], xf2*a0);
        xf1=xf2; xf2=xf3;
      }

      xma = xf[j]/a0; // 'real' sin(ph)

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

      if (fabs(xma)<level) {
        double y = 2*M_PI*period + (slope>0 ? asin(xma):(M_PI-asin(xma))) - pp; // dp
        if (!isnan(y)){ pS++; pSx+=t; pSxx+=t*t; pSxxx+=t*t*t; pSxxxx+=t*t*t*t;
                        pSy+=y; pSxy+=t*y; pSxxy+=t*t*y; }
        if (F2 && print_data==2)  fprintf(F2, "%14e %14e\n", (i+j)*dt, xp);
        if (F2 && print_data==21) fprintf(F2, "%14e %14e\n", (i+j)*dt, y);
        if (F2 && print_data==22) fprintf(F2, "%14e %14e\n", (i+j)*dt, slope*a0*level);
      }
    }

    // linear fit: dph(t) = dph0 + dph1*t
//    dp0 = (pSxy*pSx-pSxx*pSy) / (pSx*pSx-pSxx*pS);
//    dp1 = (pSx*pSy-pSxy*pS) / (pSx*pSx-pSxx*pS);
//    dp2 = 0;

    // quadratic fit: dph(t) = dph0 + dph1*t + dph2*t*t
    pdet = (pSxxxx*pSxx*pS + 2*pSxxx*pSxx*pSx - pSxx*pSxx*pSxx - pSxxx*pSxxx*pS - pSx*pSx*pSxxxx);
    dp0 =  (pSxxxx*pSxx*pSy + pSxxx*pSxx*pSxy + pSxxx*pSxxy*pSx - pSxx*pSxx*pSxxy - pSxxx*pSxxx*pSy - pSxy*pSx*pSxxxx)/pdet;
    dp1 = -(pSxxxx*pSx*pSy + pSxx*pSxx*pSxy + pSxxx*pSxxy*pS - pSxxy*pSxx*pSx - pSxxx*pSxx*pSy - pSxy*pS*pSxxxx)/pdet;
    dp2 =  (pSxxx*pSx*pSy + pSxx*pSxy*pSx + pSxxy*pSxx*pS - pSxxy*pSx*pSx - pSxx*pSxx*pSy - pS*pSxy*pSxxx)/pdet;

    // second loop - adjast amp
    for (j=0; j<np; j++){
      double t = j*dt;
      double pp = p0 + 2*M_PI*f0*t + dp0 + dp1*t + dp2*t*t;
      double xp = a0*sin(pp);

      if (fabs(xp/a0) >= level) {
        double y = (xf[j]-xp)/sin(pp); // da
        aS++; aSx+=t; aSy+=y; aSxx+=t*t; aSxy+=t*y;
      }
      if (F2 && print_data==3) fprintf(F2, "%14e %14e\n", (i+j)*dt, xp);
    }

    // only mean value, da(t) = da0
//    da0 = aSy/aS;
//    da1 = 0;

    // linear fit da(t) = da0 + da1(t)
    da0 = (aSxy*aSx-aSxx*aSy) / (aSx*aSx-aSxx*aS);
    da1 = (aSx*aSy-aSxy*aS) / (aSx*aSx-aSxx*aS);

    // da(t) = da1*t;
//    da0 = 0;
//    da1 = aSxy/aSxx;

    // calculate error and print filtered signal
    err=0; aS=0;
    for (j=0; j<np; j++){
      double t = j*dt;
      double pp = p0 + 2*M_PI*t*f0 + dp0 + dp1*t + dp2*t*t;
      double xp = (a0+da0+da1*t)*sin(pp);
      err+=(xf[j]-xp) * (xf[j]-xp); aS++;
      if (F2 && print_data==4) fprintf(F2, "%14e %14e\n", (i+j)*dt, xp);
    }
    err = sqrt(err/aS);
    if (print_data>0) fprintf(F2, "\n");

    fprintf(stderr, "%d: dph = %e + %e t + %e t^2, p0 -> %e\n", i, dp0, dp1, dp2, p0);
    fprintf(stderr, "damp = %e + %e t, a0 -> %e\n", da0, da1, a0);

    if (err > a0+da0 + 0.5*j*dt * da1) {state=0; i+=np; continue;}

    // output data for the middle of window
    fprintf(F3, "%14e %14e %14e %14e %14e\n",
      (i+0.5*j)*dt,
      p0 + dp0 + 0.5*j*dt * dp1 + 0.25*dp2*j*dt*j*dt,
      f0 + (dp1 + dp2*j*dt)/2/M_PI,
      a0 + da0 + 0.5*j*dt * da1,
      err);

    // initial values for the next step
    p0 += dp0 + j*dt * dp1  + dp2*j*dt*j*dt;
    f0 += (dp1 + 2*dp2*j*dt)/2/M_PI;
    a0 += da0 + j*dt * da1;

    i += np;
  }

  fclose(F1);
  fclose(F2);
  fclose(F3);

  return 0;
}
