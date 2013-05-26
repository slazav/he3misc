#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// number of points
#define PTS 50000

int
main(){
  int i;
  double dt=1.5e-6;
  double np = 30;   // window (number of periods)
  double level=0.5;

  double xm[PTS];

  FILE *F1 = fopen("data1.dat", "w");
  FILE *F2 = fopen("data2.dat", "w");
  FILE *F3 = fopen("data3.dat", "w");

  if (!F1 || !F2) return 1;

  // model signal
  {
    double amp0=1.2;     // signal amplitude
    double amp0n = 0.12; // noise amplitude
    double fre0=34567.0; // final frequency
    double df0=1234.0;   // frequency change
    double ta=0.1;       // amplitude relaxation time
    double tf=0.021;     // frequency relaxation time
    double ph=0.4;         // initial phase

    for (i=0; i<PTS; i++){
      double t=i*dt;
      double noise = amp0n * (double)random() / RAND_MAX - amp0n/2;
      double fre = fre0 + df0 * exp(-t/tf);
      double amp = amp0 * exp(-t/tf);
      if (i>0) ph=ph+2*M_PI*fre*dt;
      while (ph>M_PI) ph-=2*M_PI;
      xm[i] = noise + amp * sin(ph);
      fprintf(F1, "%14e %14e %14e %14e %14e\n", t, xm[i], fre, amp, amp0n);
    }
  }

  {
  // initial conditions
  double f0 = 34567.0+1000; // initial conditions
  double a0 = 1.2;
  double p0 = 0;

  int filter=1; // do freq filtering
  int print_data = 0; // print debugging information to data2.dat
                      // 0  - nothing
                      // 1  - filtered signal
                      // 2  - calculated signal before adjusting phase
                      // 21 - phase adjastment
                      // 22 - period tracking during phase adjastment
                      // 3  - calculated signal before adjusting amp
                      // 4  - final calculated signal

  // filtering
  int ipi2 = 1/f0/dt; // number of points per period
  double xf1=xm[ipi2/4], xf2=xm[ipi2/4+1]; // initial conditions for filtering oscillator
                                           // (-pi/2 phase shift)
  // phase locking
  double slope = ( (xf2>xf1 || xf2<-level) && xf2<level)? 1:-1; // initial slope
  double slope_lock; // for locking slope at short times
  double period = 0; // signal period

  for (i=0; i<PTS; i++){
    int j;

    double aS=0, aSx=0, aSy=0, aSxx=0, aSxy=0;
    double pS=0, pSx=0, pSxx=0, pSxxx=0, pSxxxx=0, pSy=0, pSxy=0, pSxxy=0;

    double pdet, dp0,dp1,dp2,da0,da1;
    period=0;
    slope_lock=0;

    // first loop - adjast phase
    for (j=i; (j-i)*dt*f0 < np; j++){
      double t = (j-i)*dt;
      double w = 2*M_PI*f0;
      double pp = p0 + w*t;
      double xp = a0*sin(pp);
      double xma;

      if (j>=PTS) break;
      xma = xm[j]/a0; // 'real' sin(ph)

      // filtering oscillator (-pi/2 phase)
      if (filter && j>1){
        double g=0.5;
        double xf3 = (xf2*(2-dt*dt*w*w) - xf1*(1-dt*w*g) -
                      xm[j]*dt*dt*w*w)/(1+dt*w*g);
        if (print_data==1) fprintf(F2, "%14e %14e\n", j*dt, xf2);
        xma = xf2/a0;
        xf1=xf2; xf2=xf3;
      }

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
        if (print_data==2) fprintf(F2, "%14e %14e\n", j*dt, xp);
        if (print_data==21) fprintf(F2, "%14e %14e\n", j*dt, y);
        if (print_data==22) fprintf(F2, "%14e %14e\n", j*dt, slope*a0*level);
      }
    }

//    dp0 = (pSxy*pSx-pSxx*pSy) / (pSx*pSx-pSxx*pS);
//    dp1 = (pSx*pSy-pSxy*pS) / (pSx*pSx-pSxx*pS);
//    dp2 = 0;

    pdet = (pSxxxx*pSxx*pS + 2*pSxxx*pSxx*pSx - pSxx*pSxx*pSxx - pSxxx*pSxxx*pS - pSx*pSx*pSxxxx);
    dp0 =  (pSxxxx*pSxx*pSy + pSxxx*pSxx*pSxy + pSxxx*pSxxy*pSx - pSxx*pSxx*pSxxy - pSxxx*pSxxx*pSy - pSxy*pSx*pSxxxx)/pdet;
    dp1 = -(pSxxxx*pSx*pSy + pSxx*pSxx*pSxy + pSxxx*pSxxy*pS - pSxxy*pSxx*pSx - pSxxx*pSxx*pSy - pSxy*pS*pSxxxx)/pdet;
    dp2 =  (pSxxx*pSx*pSy + pSxx*pSxy*pSx + pSxxy*pSxx*pS - pSxxy*pSx*pSx - pSxx*pSxx*pSy - pS*pSxy*pSxxx)/pdet;
    if (filter) dp0 -= M_PI/2;

    // second loop - adjast amp
    for (j=i; (j-i)*dt*f0 < np; j++){
      double w = 2*M_PI*f0;
      double t = (j-i)*dt;
      double pp = p0 + w*t + dp0 + dp1*t + dp2*t*t;
      double xp = a0*sin(pp);

      if (j>=PTS) break;
      while (pp>M_PI) pp-=2*M_PI;

      if (fabs(xp/a0) >= level) {
        double y = (xm[j]-xp)/sin(pp); // da
        aS++; aSx+=t; aSy+=y; aSxx+=t*t; aSxy+=t*y;
      }
      if (print_data==3) fprintf(F2, "%14e %14e\n", j*dt, xp);
    }

//    da0 = (aSxy*aSx-aSxx*aSy) / (aSx*aSx-aSxx*aS);
//    da1 = (aSx*aSy-aSxy*aS) / (aSx*aSx-aSxx*aS);

//    da0 = aSy/aS;
//    da1 = 0;

    da0 = 0;
    da1 = aSxy/aSxx;

    // print filtered signal
    for (j=i; (j-i)*dt*f0 < np; j++){
      double t = (j-i)*dt;
      double pp = p0 + 2*M_PI*t*f0 + dp0 + dp1*t + dp2*t*t;
      double xp = (a0+da0+da1*t)*sin(pp);
      if (j>=PTS) break;
      if (print_data==4) fprintf(F2, "%14e %14e\n", j*dt, xp);
    }
    if (print_data>0) fprintf(F2, "\n");

    fprintf(F3, "%14e %14e %14e %14e\n", i*dt, p0+dp0, f0+dp1/2/M_PI , a0+da0);

    // initial values for next step
    p0 += dp0 + (j-i)*dt * dp1  + dp2*(j-i)*dt*(j-i)*dt;
    f0 += (dp1 + dp2*(j-i)*dt)/2/M_PI;
    a0 += da0 + (j-i)*dt * da1;

    fprintf(stderr, "dph = %e + %e t + %e t^2, p0 -> %e\n", dp0, dp1, dp2, p0);
    fprintf(stderr, "damp = %e + %e t, a0 -> %e\n", da0, da1, a0);

    if (j<=i+1) break;
    i = j-1;
  }
  }

  fclose(F1);
  fclose(F2);
  fclose(F3);

  return 0;
}
