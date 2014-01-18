#include <stdio.h>
#include <stdlib.h>

FILE *pnm_file;
int pnm_cols=7;
int pnm_space=5;
int pnm_lcount=0;
int pnm_npts=0;

struct arrays_t{
  double usol;
  double xsol;
};

extern struct arrays_t arrays_;

void
pnm_open_(int *npts){
  pnm_file = fopen("result.pnm", "w");
  pnm_npts = *npts;
  fprintf(pnm_file, "P6\n");
  fprintf(pnm_file, "%5d %5d\n", (pnm_npts+pnm_space)*pnm_cols,0);
  fprintf(pnm_file, "255\n");
}

void pnm_color(double x, double xmin, double xmax){

  unsigned char r,g,b;
  if (x>xmax) x=xmax;
  if (x<xmin) x=xmin;
  r = ((x-xmin)/(xmax-xmin)*255.0);
  g = ((x-xmin)/(xmax-xmin)*255.0);
  b = ((x-xmin)/(xmax-xmin)*255.0);
  fprintf(pnm_file, "%c%c%c", r,g,b);
}

void
pnm_write_(){
  int i,j;
  double xmin[7] = {-1,-1,-1, -1,-1,-1, -0.28};
  double xmax[7] = { 1, 1, 1,  1, 1, 1, -0.24};
  for (j=0; j<pnm_cols; j++){
    for (i=0; i<pnm_npts; i++){
      double val = (&arrays_.usol)[7*i + j];
      if (j==6) val=cos(val);
      pnm_color(val, xmin[j], xmax[j]);
    }
    for (i=0; i<pnm_space; i++){
      fprintf(pnm_file, "%c%c%c", 0,0,0);
    }
  }
  pnm_lcount++;
}

void
pnm_close_(){
  fseek(pnm_file, 9, SEEK_SET);
  fprintf(pnm_file, "%5d", pnm_lcount);
  fclose(pnm_file);
}
