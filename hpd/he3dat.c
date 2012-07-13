#include "he3dat.h"
#include <stdio.h>
#include <stdlib.h>

void usage(){
  printf(
    "he3dat -- print He3 parameters\n"
    "Usage: he3dat <command> <parameters>\n"
    "\n"
    "Commands:\n"
    "  tc <pressure, bar>\n"
    "  tc <p1> <p2> <step> --\n"
  );
  exit(1);
}

void error(const char * err){
  printf("Error: %s\n", err);
  exit(1);
}

int
main(int argc, char *argv[]){

  if (argc < 2) usage();

  /* Tc */
  if (strcasecmp(argv[1], "tc")==0){
    if (argc == 3){
      double p = atof(argv[2]);
      printf("%f\n", he3_tc_(&p));
    }
    else if (argc == 5){
      double p1 = atof(argv[2]);
      double p2 = atof(argv[3]);
      double st = atof(argv[4]);
      double p;
      if ( p2 <= p1 || st <= 0 )
        error("Wrong pressure range");
      for (p=p1; p<=p2; p+=st){
        printf("%f %f\n", p, he3_tc_(&p));
      }
    }
    else error("Unknown command");
    exit(0);
  }

  /* Tab */
  if (strcasecmp(argv[1], "tab")==0){
    if (argc == 3){
      double p = atof(argv[2]);
      printf("%f\n", he3_tab_(&p));
    }
    else if (argc == 5){
      double p1 = atof(argv[2]);
      double p2 = atof(argv[3]);
      double st = atof(argv[4]);
      double p;
      if ( p2 <= p1 || st <= 0 )
        error("Wrong pressure range");
      for (p=p1; p<=p2; p+=st){
        printf("%f %f\n", p, he3_tab_(&p));
      }
    }
    else error("Wrong parameters");
    exit(0);
  }

  /* Pmelt */
  if (strcasecmp(argv[1], "pmelt")==0){
    if (argc == 3){
      double t = atof(argv[2]);
      printf("%f\n", he3_pmelt_(&t));
    }
    else if (argc == 5){
      double t1 = atof(argv[2]);
      double t2 = atof(argv[3]);
      double st = atof(argv[4]);
      double t;
      if ( t2 <= t1 || st <= 0 )
        error("Wrong temperature range");
      for (t=t1; t<=t2; t+=st){
        printf("%f %f\n", t, he3_pmelt_(&t));
      }
    }
    else error("Wrong parameters");
    exit(0);
  }


}