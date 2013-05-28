#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// number of points
#define PTS 50000

int
main(){
  int i;
  double dt=1.5e-6;   // time step
  double dx=1e-4;
  FILE *F1, *F2;

  F1 = fopen("data.dat", "w");      // signal itself
  F2 = fopen("data_test.dat", "w"); // all parameters
  if (!F1 || !F2) return 1;

  {
    double amp0=1.2;     // signal amplitude
    double amp0n = 0.22; // noise amplitude
    double fre0=34567.0; // final frequency
    double df0=1234.0;   // frequency change
    double ta=0.012;     // amplitude relaxation time
    double tf=0.021;     // frequency relaxation time
    double ph=0.4;       // initial phase
    double t0=0.0005;    // start of signal
    double base=0;

    fprintf(F1, "+1,+0,+%d,+1,%+11E,%+11E,+0,%+14E,+0.0,+0\n\n",
                 PTS, dt,-t0,dx);

    fprintf(F2, "#%13s %14s %14s %14s %14s\n",
      "t", "x", "fre", "amp", "noise");
    for (i=0; i<PTS; i++){
      double t=i*dt;
      double noise = amp0n * (double)random() / RAND_MAX - amp0n/2;
      double fre = fre0 + df0 * exp(-t/tf);
      double amp = t<t0? 0 : amp0 * exp(-t/ta);
      double x;
      if (i>0) ph=ph+2*M_PI*fre*dt;
      while (ph>M_PI) ph-=2*M_PI;
      x = noise + amp * sin(ph) + base;
      fprintf(F1, "%d\n", (int)(x/dx));
      fprintf(F2, "%14e %14e %14e %14e %14e\n", t-t0, x, fre, amp, amp0n/2);
    }
  }

  fclose(F1);
  fclose(F2);

  return 0;
}
