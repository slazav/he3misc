
MODULE general

  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(8)

  INTEGER, PARAMETER :: maxnpt = 10000

  !these are initializations, values are reset in mainlib.f90
  INTEGER, SAVE :: nmax=250
  INTEGER, SAVE :: nrmax=250
  INTEGER, SAVE :: nfmax=250
  INTEGER, SAVE :: nmaxt=250

  REAL (KIND=dp), PARAMETER :: pi=3.14159265358979_dp

END MODULE
