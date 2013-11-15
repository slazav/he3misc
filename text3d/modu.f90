MODULE modu

  USE general
  USE text

  !! parameters used to calculate texture

  REAL (KIND=dp), SAVE :: t,p,h,r,lo, zrange
  REAL (KIND=dp), SAVE :: chia, vd, xir, de, dar, lsg
  REAL (KIND=dp), SAVE :: dx, stilt
  REAL (KIND=dp), SAVE :: nub,nu0

  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: apsi
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evr,evf,evz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elr,elf,elz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: ew


  contains

  subroutine set_text_pars(tt,pp,hh)
    REAL (KIND=dp) :: tt, pp, hh
    chia = fchia(tt,pp)
    vd   = fvd(tt,pp)
    xir  = fxih(tt,pp,hh)/r
    de   = fdelta(tt,pp)
    dar  = fdar(tt,pp,r)
    lsg  = 3._dp ! see Fig. 1 in Erkki's paper
  end subroutine

END MODULE


