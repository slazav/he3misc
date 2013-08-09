extern struct{
   double he3_gyro_, he3_amass_, he3_mmass_;
   double const_na_, const_kb_, const_r_, const_h_, const_hbar_;
   double const_pi_;
   double he3_pcr_,
          he3_tcr_,
          he3_pabn_,
          he3_tabn_,
          he3_psmin_,
          he3_tsmin_,
          he3_pa_,
          he3_ta_,
          he3_pb_,
          he3_tb_,
          he3_pneel_,
          he3_tneel_;
   /* order in important - see common block in he3.fh */
} he3_const_;

double he3_pvap_(double *T);   /* Vapor pressure [bar] vs T [K] */
double he3_pmelt_(double *T);  /* Melting pressure [bars] vs T [K] */
double he3_tc_(double *P);     /* T_c [mK] vs P [bar] */
double he3_tab_(double *P);    /* T_ab [mK] vs P [bar] */

double he3_vm_(double *P);     /* Molar Volume [cm**3/mole] vs P [bar] */
double he3_gammaf_(double *P); /* R-Gas constant GAMMA=C/RT [1/(K*mol)] vs P [bar] */
double he3_c1_(double *P);     /* First sound velosity [cm/s] vs P [bar] */
double he3_tmag_(double *P);   /* Magnetic temperature [K] vs P [bar] */

double he3_rho_(double *P);    /* Density [g/cm^3] vs P [bar] */
double he3_2n0_(double *P);    /* 2N0 vs P [bar] */
double he3_mm_(double *P);     /* Effective mass / atom mass vs P [bar] */
double he3_meff_(double *P);   /* Effective mass [g] vs P [bar] */
double he3_pf_(double *P);     /* Fermi momentum [sgs] vs P [bar] */
double he3_vf_(double *P);     /* Fermi velocity [cm/s] vs P [bar] */
double he3_chi0_(double *P);   /* Susceptibility vs P [bar]  */
double he3_f0a_(double *P);    /* F0a vs P [bar] */
double he3_f0s_(double *P);    /* F0s vs P [bar] */
double he3_f1a_(double *P);    /* F1a vs P [bar] */
double he3_f1s_(double *P);    /* F1s vs P [bar] */
double he3_z0_(double *P);     /* Z0 vs P [bar]  */
double he3_a_(double *P);      /* average atomic spacing, angstr. */
double he3_gdk_(double *P);    /* average dipolar coupling energy, K */
double he3_tfeff_(double *P);  /* effective fermi temperature, K */

double he3_flegg_(double *P, double *ttc); /* Legget freq^2, [Hz^2] vs P [bar], T/Tc */

double he3_swvel_(double *P, double *ttc);     /* Osheroff's spin wave vel. [cm/s] vs P [bar], T [mK] */
double he3_swvel_par_(double *P, double *ttc); /* Perp Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */
double he3_swvel_per_(double *P, double *ttc); /* Parallel Fomin spin wave vel. [cm/c] vs P [bar], T [mK] */

double he3_ds_exp_(double *P, double *ttc); /* spin diffusion coeff. in superfluid He3 (measured) */
double he3_dn_exp_(double *P, double *ttc); /* spin diffusion in normal He3 */
double he3_d_exp_(double *P, double *ttc);  /* combined normal + superfluid spin diffusion */

double he3_susept_(double *P, double *ttc); /* Suseptibility [sgs] vs P [bar], T [mK] */

double he3_tau_r_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */
double he3_tau_f_(double *ttc);     /* Leggett-Takagi tau_r [s] vs T/Tc, 20bar */

double he3_yosida_(double *ttc);    /* Yosida function vs T/Tc */