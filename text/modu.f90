MODULE modu

  USE general

  REAL (KIND=dp), SAVE :: t,p,h,r,lo

  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evr,evf,evz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elr,elf,elz
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: ew

END MODULE


