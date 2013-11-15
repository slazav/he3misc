
MODULE general

  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(8)
  INTEGER, PARAMETER :: qp=SELECTED_REAL_KIND(17)
  INTEGER, PARAMETER :: maxnpt = 500000

  
  !these are initializations, values are reset in mainlib.f90
  INTEGER, SAVE :: nmax=250
  INTEGER, SAVE :: nrmax=250
  INTEGER, SAVE :: nfmax=250
  INTEGER, SAVE :: nmaxt=250
  INTEGER, SAVE :: nzmax=250
  INTEGER, SAVE :: nmaxp=250

  INTEGER, SAVE :: interpolation_error=0

  REAL (KIND=dp), PARAMETER :: pi=3.14159265358979_dp

END MODULE
