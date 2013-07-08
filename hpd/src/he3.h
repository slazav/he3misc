extern struct{
   double he3_pabn_, he3_tabn_, he3_pa_;
   double he3_amass_, he3_gyro_;
   double ana, r, hc, akb, pi;
   /* order in important - see common block in he3.fh */
} he3_const_;

double he3_pmelt_(double *T);  /* Melting pressure [bars] vs T [mK] */
double he3_Tc_(double *P);     /* T_c [mK] vs P [bar] */
double he3_tab_(double *P);    /* T_ab [mK] vs P [bar] */

double he3_vm_(double *P);     /* Molar Volume [cm**3/mole] vs P [bar] */
double He3_Meff_(double *P);   /* Effective mass [g] vs P [bar] */
double He3_Pf_(double *P);     /* Fermi momentum [sgs] vs P [bar] */
double He3_Vf_(double *P);     /* Fermi velocity [cm/s] vs P [bar] */

double He3_gammaf_(double *P); /* R-Gas constant GAMMA=C/RT [1/(K*mol)] vs P [bar] */
/*
double He3_dnde      ! Density of state

double He3_Flegg     ! Legget freq^2, [Hz^2] vs P [bar], T/Tc

double He3_swvel     ! Osheroff's spin wave vel. [cm/s] vs P [bar], T [mK]
double he3_swvel_par ! Perp Fomin spin wave vel. [cm/c] vs P [bar], T [mK]
double he3_swvel_per ! Parallel Fomin spin wave vel. [cm/c] vs P [bar], T [mK]

double He3_Ds_exp    ! spin diffusion coeff. in superfluid He3 (measured)
double He3_Dn_exp    ! spin diffusion in normal He3
double He3_D_exp     ! combined normal + superfluid spin diffusion

double He3_susept    ! Suseptibility [sgs] vs P [bar], T [mK]

double He3_tau_r     ! Leggett-Takagi tau_r [s] vs T/Tc, 20bar
double He3_tau_f     ! Leggett-Takagi tau_r [s] vs T/Tc, 20bar

double He3_yosida    ! Yosida function vs T/Tc
double He3_Z0        ! Z0 vs P [bar]
double He3_F0a       ! F0A vs P [bar]
*/