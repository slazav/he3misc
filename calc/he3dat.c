#include "he3dat.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>

/* max number of table columns >= 4 */
#define MAXCOL 10

/* Print usage information and exit */
void
usage(){
  printf(
    "he3dat -- print He3 parameters\n"
    "Usage: he3dat <options> <command>\n"
    "\n"
    "Options:\n"
    "  --help, -h     -- print usage information and exit\n"
    "  --verb, -v     -- print additional information\n"
    "  --temp, -t <value>  -- set temperature value or range [mK]\n"
    "  --pres, -p <value>  -- set pressure value or range [bar]\n"
    "  --ttc,  -c <value>  -- use T/Tc instead of T (pressure must be set)\n"
    "\n"
    "Commands (set -p ot -t options ):\n"
    "  Tabn, Pabn    -- A-B-Normal critical point\n"
    "  Pa            -- A-Solid-Normal critical point\n"
    "  Pmelt(T)      -- Melting pressure\n"
    "  Tc(P), Tab(P), Tab/Tc(P)  -- Tc, Tab, Tab/Tc\n"
    "  Vm(P)         -- Molar Volume\n"
    "  Meff(P)       -- Effective mass\n"
    "  Vf(P), Pf(P)  -- Fermi velosity and momentum\n"
    "\n"
    "Examples:\n"
    "$ he3dat Vm -p 22.3     -- print molar volume for P=22.3 bar\n"
    "$ he3dat Tab/Tc -p 22:29:1 -- print Tab/Tc table for 22..29 bar\n"
    "$ he3dat Pmelt -t 0.5:1.5 -- print Melting pressure (default step 0.1 mK)\n"
  );
  exit(1);
}

/* Print short error message and exit */
void
error(const char * err){
  printf("Error: %s\n", err);
  exit(1);
}

/* Parse range setting in form <v1>:<v2>:<step>
   Return value: 0 for single value mode; 1 for table mode. */
int
parse_range(const char * str,
      double * v1, double * v2, double *st){
  const char *del=":";
  char *tok;

  tok = strtok((char*)str,  del);
  if (tok) *v1 = atof(tok);
  else error("bad parameter");

  tok = strtok(NULL, del);
  if (tok) *v2 = atof(tok);
  else return 0;

  tok = strtok(NULL, del);
  *st = tok? atof(tok):0.1;

  if (*st <= 0 || *v2 < *v1)
    error("bad table parameters");

  return 1;
}

/* Some commands needs no parameters, T or P or T and P.
  This function check this... */
void check_pars(int ptab, int pneed, int ttab, int tneed){
  if (ptab<0 && pneed) error("pressure must be set; use -p option.");
  if (ttab<0 && tneed) error("temperature must be set; use -t option.");
  if (ptab>=0 && !pneed) error("no need for pressure option.");
  if (ttab>=0 && !tneed) error("no need for temperature option");
}

/* MAIN */
int
main(int argc, char *argv[]){

  /* Parameters: -1: no par, 0: normal mode, 1: table mode.
     p - pressure, t - temp, c - t/tc */
  int ttab=-1, ptab=-1, ctab=-1, tctab;
  double t1,t2,ts, p1,p2,ps, c1,c2,cs;
  int tn=1, pn=1; /* number of table rows */
  int ti, pi;
  double *table[MAXCOL];     /* data table */
  const char *thead[MAXCOL]; /* column heads */
  int pars; /* number of calculated columns 0..MAXCOL-3*/

  int v=1; /* verbosity level */
  int c, i;
  const char *cmd;

  /* Parse command line options */
  if (argc==1) usage();
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"quiet", no_argument, 0,  'q' },
      {"help",  no_argument, 0,  'h' },
      {"temp",  required_argument, 0,  't' },
      {"pres",  required_argument, 0,  'p' },
      {"ttc",   required_argument, 0,  'c' },
      {0,0,0,0} };
    c = getopt_long(argc, argv, "qht:p:c:",
        long_options, &option_index);
    if (c == -1)  break;
    switch (c) {
    case 'c': ctab = parse_range(optarg, &c1,&c2,&cs); break;
    case 't': ttab = parse_range(optarg, &t1,&t2,&ts); break;
    case 'p': ptab = parse_range(optarg, &p1,&p2,&ps); break;
    case 'q': v=0; break;
    case 'h': usage();
    default: exit(1);
    }
  }

  /* we can't use both T and T/Tc settings */
  if (ttab>=0 && ctab>=0) error("can't use both T and T/Tc settings");
  tctab = ttab>=0 ? ttab : ctab;

  /* we can't use T/Tc if we don't know pressure */
  if (ctab>=0 && ptab<0) error("can't use T/Tc without P setting");

  /* create table and fill P, T, TTC columns */
  if (ptab>0) pn = ceil((p2-p1)/ps) + 1;
  if (ttab>0) tn = ceil((t2-t1)/ts) + 1;
  if (ctab>0) tn = ceil((c2-c1)/cs) + 1;
  for (i=0; i<MAXCOL; i++)
    table[i] = (double *) malloc(pn*tn * sizeof(double));
  for (pi = 0; pi< pn; pi++){
    for (ti = 0; ti< tn; ti++){
      table[0][pi*tn+ti] = -1; /* P */
      table[1][pi*tn+ti] = -1; /* T */
      table[2][pi*tn+ti] = -1; /* TTC */

      if (ptab==0) table[0][pi*tn+ti] = p1;
      if (ttab==0) table[1][pi*tn+ti] = t1;
      if (ctab==0) table[2][pi*tn+ti] = c1;

      if (ptab==1) table[0][pi*tn+ti] = p1 + ps * pi;
      if (ttab==1) table[1][pi*tn+ti] = t1 + ts * ti;
      if (ctab==1) table[2][pi*tn+ti] = c1 + cs * ti;

      /* calculate T/Tc from T or T from T/Tc */
      if (ctab>=0 && ptab>=0) table[1][pi*tn+ti] =
           table[2][pi*tn+ti] * he3_tc_(&table[0][pi*tn+ti]);
      if (ttab>=0 && ptab>=0) table[2][pi*tn+ti] =
           table[1][pi*tn+ti] / he3_tc_(&table[0][pi*tn+ti]);
    }
  }


  /* we need 1..MAXCOL-3 command parameter */
  pars=argc-optind;
  if (pars < 1 ) error("single command needed");
  if (pars > MAXCOL-3 ) error("too many commands");

  /* fill the rest of the table */
  for (i=0; i<pars; i++){
    cmd=argv[optind+i];

    /* parse commands */
    if (strcasecmp(cmd, "Pa")==0){
      check_pars(ptab, 0, tctab, 0);
      thead[i+3] = "Pa, bar";
      table[i+3][0] = he3consts_.he3_pa;
    }
    else if (strcasecmp(cmd, "Pabn")==0){
      check_pars(ptab, 0, tctab, 0);
      thead[i+3] = "Pabn, bar";
      table[i+3][0] = he3consts_.he3_pabn;
    }
    else if (strcasecmp(cmd, "Tabn")==0){
      check_pars(ptab, 0, tctab, 0);
      thead[i+3] = "Tabn, mK";
      table[i+3][0] = he3consts_.he3_tabn;
    }
    else if (strcasecmp(cmd, "Pmelt")==0){
      check_pars(ptab, 0, tctab, 1);
      if (v) he3_pmelt_comm_();
      thead[i+3] = "Pmelt, bar";
      for (ti=0; ti<tn; ti++)
        table[i+3][ti] = he3_pmelt_(&table[1][ti]);
    }
    else if (strcasecmp(cmd, "Tc")==0){
      check_pars(ptab, 1, tctab, 0);
      if (v) he3_tc_comm_();
      thead[i+3] = "Tc, mK";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_tc_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Tab")==0){
      check_pars(ptab, 1, tctab, 0);
      if (v) he3_tab_comm_();
      thead[i+3] = "Tab, mK";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_tab_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Tab/Tc")==0){
      check_pars(ptab, 1, tctab, 0);
      thead[i+3] = "Tab/Tc";
      for (pi=0; pi<pn; pi++) table[i+3][pi] =
              he3_tab_(&table[0][pi])/he3_tc_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Vm")==0){
      check_pars(ptab, 1, tctab, 0);
      if (v) he3_vm_comm_();
      thead[i+3] = "Vm, cm3/mole";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_vm_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Meff")==0){
      check_pars(ptab, 1, tctab, 0);
      thead[i+3] = "Meff, g";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_meff_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Pf")==0){
      check_pars(ptab, 1, tctab, 0);
      thead[i+3] = "Pf, g cm/s";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_pf_(&table[0][pi]);
    }
    else if (strcasecmp(cmd, "Vf")==0){
      check_pars(ptab, 1, tctab, 0);
      thead[i+3] = "Vf, cm/s";
      for (pi=0; pi<pn; pi++)
        table[i+3][pi] = he3_vf_(&table[0][pi]);
    }

    else error("unknown command");
  }

  /* print table title */
  if (v){
    printf("#\n");
    if (ptab==0)
      printf("# Pressure: %.4f bar\n", p1);
    else if (ptab==1)
      printf("# Pressure: %.4f .. %.4f bar, step %.4f bar\n", p1, p2, ps);

    if (ttab==0)
      printf("# Temperature: %.4f mK\n", t1);
    else if (ttab==1)
      printf("# Temperature: %.4f .. %.4f mK, step %.4f mK\n", t1, t2, ts);

    if (ctab==0)
      printf("# T/Tc: %.4f\n", c1);
    else if (ctab==1)
      printf("# T/Tc: %.4f .. %.4f, step %.4f\n", c1, c2, cs);
  }

  /* print column titles */
  if (v){
    printf("#\n");
    printf("# ");
    /* print pressure column if we have pressure range */
    if (ptab==1) printf("%-7s ", "P, bar");
    /* print temp column if we have T or T/Tc range */
    if (tctab==1) printf("%-7s ", "T, mK");
    /* print T/Tc column if we have T or T/Tc range and P is set*/
    if (ptab>=0 && tctab==1) printf("%-7s ", "T/Tc");

    /* print other titles */
    for (i=0; i<pars; i++) printf("%-12s ", thead[i+3]);
    printf("\n");
  }

  /* print table! */
  for (pi = 0; pi< pn; pi++){
    for (ti = 0; ti< tn; ti++){
      double P = table[0][pi*tn+ti];
      double T = table[1][pi*tn+ti];
      double C = table[2][pi*tn+ti];
      if (P<0) P=nan("");
      if (T<0) T=nan("");
      if (C<0) C=nan("");

      if (ptab==1) printf("%7.4g ", P);
      if (tctab==1) printf("%7.4g ", T);
      if (ptab>=0 && tctab==1) printf("%7.4g ", C);

      printf("  "); /* to compensate leading "# " in the header */

      for (i=0; i<pars; i++)
        printf("%-12.6g ", table[i+3][pi*tn+ti]);

      printf("\n");
    }
    if (pi!=pn-1 && tn > 1) printf("\n");
  }

  exit(0);
}