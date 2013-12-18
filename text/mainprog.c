#include <stdio.h>
#include "texture.h"

int
main(){
  FILE *F;
  int npt=200, ns=500;
  int i, initype,msglev=1;
  struct text_struct pars;

  double textpar[10], specpar[2];
  double textur[3][npt+1], spec[2][ns], apsi[npt+1];

  F = fopen("initials.dat", "r");
  fscanf(F, "%lf %*[^\n]", &pars.ttc); // t
  fscanf(F, "%lf %*[^\n]", &pars.p); // p
  fscanf(F, "%lf %*[^\n]", &pars.f0); // nu0
  fscanf(F, "%lf %*[^\n]", &pars.r); // r
  fscanf(F, "%lf %*[^\n]", specpar+0); // gamma
  fscanf(F, "%lf %*[^\n]", specpar+1); // fac
  fscanf(F, "%lf %*[^\n]", &pars.omega); // omega
  fscanf(F, "%lf %*[^\n]", &pars.omega_v); // ov
  fscanf(F, "%lf %*[^\n]", &pars.lo); // lo
  fscanf(F, "%i %*[^\n]",  &initype);  //
  fscanf(F, "%lf %*[^\n]", &pars.lhv); // lhv
  fscanf(F, "%lf %*[^\n]", &pars.chi); // chi
  fscanf(F, "%lf %*[^\n]", &pars.nuB); // nuB
  fclose(F);

  for (i=0; i<=npt; i++) apsi[i] = 0.0;
  calctexture_(&npt,&pars, &ns,&specpar,&initype,&textur,&spec,&msglev,&apsi);

  F = fopen("texture.dat", "w");
  fprintf(F, "#r, alpha, beta\n");

  // save texture
  for (i=0; i<=npt; i++)
    fprintf(F, "%le %le %le\n", textur[0][i], textur[1][i], textur[2][i]);
  fclose(F);

  // Save the NMR spectrum
  F = fopen("spec.dat", "w");
  for (i=0; i<ns; i++)
    fprintf(F, "%le %le\n", spec[0][i], spec[1][i]);
  fclose(F);
}
