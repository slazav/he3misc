MODULE modu

  USE general

  !! parameters used to calculate texture

  REAL (KIND=dp), SAVE :: t,p,h,r,lo
  REAL (KIND=dp), SAVE :: chi, flhvfix, chia, vd, xir, de, dar, lsg
  REAL (KIND=dp), SAVE :: dx
  REAL (KIND=dp), SAVE :: nub,nu0

  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: apsi
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evr,evf,evz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elr,elf,elz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: ew


  contains

  subroutine set_text_pars(tt,pp,hh)
    REAL (KIND=dp) :: tt, pp, hh
    include '../lib/he3.f90h'
    chia = chi/(he3_text_a(t,p))
    vd   = he3_text_vd(tt,pp)
    xir  = he3_text_xih(tt,pp,hh)/r
    de   = he3_text_delta(tt,pp)
    dar  = he3_text_d(t,p) / (he3_text_a(t,p)*r)
    lsg  = 3._dp ! see Fig. 1 in Erkki's paper
  end subroutine

END MODULE


