extern double he3_pmelt_(double *t); /* Melting Pressure */
extern double he3_pmelt_comm_();
extern double he3_tc_(double *p);    /* Tc */
extern double he3_tc_comm_();
extern double he3_tab_(double *p);   /* Tab */
extern double he3_tab_comm_();
extern double he3_vm_(double *p);    /* molar volume */
extern double he3_vm_comm_();

extern double he3_meff_(double *p);  /* effective mass */
extern double he3_pf_(double *p);    /* Fermi momentum */
extern double he3_vf_(double *p);    /* Fermi velosity */

extern struct he3const_t{
  double he3_pabn;
  double he3_tabn;
  double he3_pa;
  double ana, hc, r, akb, gam, am3;
} he3consts_;
