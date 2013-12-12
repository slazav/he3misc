#include <stdio.h>

int
main(){
  FILE *F;
  int npt=200, ns=500;
  int i, initype,msglev=1;

  double textpar[10], specpar[2];
  double textur[3][npt+1], spec[2][ns+1], apsi[npt+1];

  F = fopen("initials.dat", "r");
  fscanf(F, "%lf %*[^\n]", textpar+0); // t
  fscanf(F, "%lf %*[^\n]", textpar+1); // p
  fscanf(F, "%lf %*[^\n]", textpar+2); // nu0
  fscanf(F, "%lf %*[^\n]", textpar+3); // r
  fscanf(F, "%lf %*[^\n]", specpar+0); // gamma
  fscanf(F, "%lf %*[^\n]", specpar+1); // fac
  fscanf(F, "%lf %*[^\n]", textpar+4); // omega
  fscanf(F, "%lf %*[^\n]", textpar+5); // ov
  fscanf(F, "%lf %*[^\n]", textpar+6); // lo
  fscanf(F, "%i %*[^\n]",  &initype);  //
  fscanf(F, "%lf %*[^\n]", textpar+7); // lhv
  fscanf(F, "%lf %*[^\n]", textpar+8); // chi
  fscanf(F, "%lf %*[^\n]", textpar+9); // nuB
  fclose(F);

  for (i=0; i<=npt; i++) apsi[i] = 0.0;
  calctexture_(&npt,&textpar,&ns,&specpar,&initype,&textur,&spec,&msglev,&apsi);

  F = fopen("texture.dat", "w");
  fprintf(F, "#r, alpha, beta\n");

  // save texture
  for (i=0; i<=npt; i++)
    fprintf(F, "%le %le %le\n", textur[0][i], textur[1][i], textur[2][i]);
  fclose(F);

  // Save the NMR spectrum
  F = fopen("spec.dat", "w");
  for (i=1; i<=ns; i++)
    fprintf(F, "%le %le\n", spec[0][i], spec[1][i]);
  fclose(F);
}
