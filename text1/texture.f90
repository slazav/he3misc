program texture
  implicit none

  write(*,*) 'test'

end

!!! Solve variation problem with fixed step H, and number of points
!  U     -- on input: initial values, on output: result
!  NU    -- U-function dimensions
subroutine VAR2D(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,H)
  integer NX,NY,NU, IX,IY,IU
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX,NY,NU), H

  !! fix initial values to fit limits and boundary conditions

  !! 
  contains
    function INT_IJ(I, J) result (IIJ)
  end
end

!!! User-defined functional for minimization
!
! 2D-texture energies in rotating 3He-B
! as a function of alpha and beta angles of n-vector
! (Gradient energy +  Orientation energy)
function VAR_FNC(r,z,U,Ur,Uz, NU) result(F)
  integer NU
  double precision r,z,U(NU),Ur(NU),Uz(NU)

  double precision CA,CB, SA,SB ! \cos\alpha, \sin\beta, ...
  double precision Nr,Nf,Nz, Nrr,Nfr,Nzr, Nrz,Nfz,Nzz
  double precision NxZ2, divN, NrotN, Fh, Fg, F

  CA = cos(U(1))
  SA = sin(U(1))
  CB = cos(U(2))
  SB = sin(U(2))

  !! Calculate n-vector and its derivatives from alpha and beta
  Nr = SB*CA ! n_r = \sin\beta \cos\alpha
  Nf = SB*SA ! n_f = \sin\beta \sin\alpha
  Nz = CB    ! n_z = \cos\beta

  Nrr = CA*CB*Ur(2)-SA*SB*Ur(1)  ! dn_r/dr
  Nrz = CA*CB*Uz(2)-SA*SB*Uz(1)  ! dn_r/dz
  Nfr = CB*SA*Ur(2)+CA*SB*Ur(1)  ! dn_f/dr
  Nfz = CB*SA*Uz(2)+CA*SB*Uz(1)  ! dn_f/dr
  Nzr = -SB*Ur(2) ! dn_z/dr = -\sin\beta * d\beta/dr
  Nzz = -SB*Uz(2) ! dn_z/dz = -\sin\beta * d\beta/dz

  !! Orientation energy:
  NxZ2 = Nf**2 + Nr**2
  Fh =  (1 - Lambda - KappaH) * NxZ2 &
         + 5D0/8D0 * Lambda * NxZ2**2

  !! Gradient energy
  divN  = Nrr + Nr/R + Nzz
  NrotN = - Nr*Nfz + Nf*(Nrz-Nzr) + Nz*(Nf/R + Nfr)
  Fg = XiH**2 * (DivN**2 - 1D0/16D0*(sqrt(3D0)*DivN + sqrt(5D0)*NRotN )**2)

  !! Full energy
  F  = Fh + Fg
end

!!! User defined texture limits and boundary conditions
subroutine VAR_LIM(r,z,Umin,Umax, NU)
  integer NU
  double precision r,z,Umin(NU),Umax(NU)

  !! limits for \alpha and \beta
  Umin(1) = -PI
  Umax(1) = PI
  Umin(2) = 0D0
  Umax(2) = 2D0*PI

  !! boundary conditions
  if (r.le.0) then
    Umin(2) = 0D0
    Umax(2) = Umin(2)
  endif
  if (r.ge.1) then
    Umin(1) = PI/3D0
    Umin(2) = Umin(1)
    Umin(2) = acos(1D0/sqrt(5D0))
    Umax(2) = Umin(2)
  endif
end

