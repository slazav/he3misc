MODULE modu

  USE general

  !! parameters used to calculate texture

  REAL (KIND=dp), SAVE :: t,p,h,r,lo
  REAL (KIND=dp), SAVE :: chi, lhv, chia, vd, xir, de, dar, lsg
  REAL (KIND=dp), SAVE :: dx
  REAL (KIND=dp), SAVE :: nub,nu0

  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: apsi
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evr,evf,evz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elr,elf,elz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: ew

END MODULE


