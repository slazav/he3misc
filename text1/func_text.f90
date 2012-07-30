!!! User-defined functional for minimization
!
! 2D-texture energies in rotating 3He-B
! as a function of alpha and beta angles of n-vector
! (Gradient energy +  Orientation energy)
subroutine VAR_FNC(r,z,U,Ur,Uz, NU, F)
  implicit none
  integer NU
  double precision r,z,U(NU),Ur(NU),Uz(NU),F

  common /texture_pars/ XiH, Lambda, KappaH
  double precision XiH, Lambda, KappaH

  double precision CA,CB, SA,SB ! \cos\alpha, \sin\beta, ...
  double precision CTA,CTB ! \ctan\alpha, \ctan\beta
  double precision Nr,Nf,Nz, Nrr,Nfr,Nzr, Nrz,Nfz,Nzz, Nxx,Nyy
  double precision NxZ2, divN, gN2, NrotN, Fh, Fg

  CA = dcos(U(1))
  SA = dsin(U(1))
  CB = dcos(U(2))
  SB = dsin(U(2))

  !! Calculate n-vector and its derivatives from alpha and beta
  Nr = SB*CA ! n_r = \sin\beta \cos\alpha
  Nf = SB*SA ! n_f = \sin\beta \sin\alpha
  Nz = CB    ! n_z = \cos\beta

  !! Orientation energy:
  Fh =  (1D0 - Lambda - KappaH) * SB**2 &
         + 5D0/8D0 * Lambda * SB**4

  Nrr = CA*CB*Ur(2)-SA*SB*Ur(1)  ! dn_r/dr
  Nrz = CA*CB*Uz(2)-SA*SB*Uz(1)  ! dn_r/dz
  Nfr = CB*SA*Ur(2)+CA*SB*Ur(1)  ! dn_f/dr
  Nfz = CB*SA*Uz(2)+CA*SB*Uz(1)  ! dn_f/dr
  Nzr = -SB*Ur(2) ! dn_z/dr = -\sin\beta * d\beta/dr
  Nzz = -SB*Uz(2) ! dn_z/dz = -\sin\beta * d\beta/dz

  ! Cf = \cos\phi, Sf=\sin\phi
  ! Nx = Nr Cf - Nf Sf
  ! Ny = Nr Sf + Nf Cf
  ! r = sqrt(x^2 + y^2), \phi = atan(y/x)
  ! dr/dx = Cf,    dr/dy=Sf
  ! df/dx = -Sf/r, df/dy = Cf/r

  ! Nxx = Nrr Cf^2 - Nfr Cf Sf - (Nrf Cf Sf - Nr Sf^2 - Nff Sf^2 - Nf Cf Sf)/r
  ! Nyy = Nrr Sf^2 + Nfr Cf Sf + (Nrf Cf Sf + Nr Cf^2 + Nff Cf^2 - Nf Cf Sf)/r

  ! divN = Nrr + Nr/r + Nff/r + Nzz

  ! in our case \phi=0, Nff=Nrf=0
  ! Nxx = Nrr, Nyy = Nr/r

  !! Gradient energy
  divN  = Nrr + Nr/r + Nzz
  gN2  = Nrr**2 + (Nr/r)**2 + Nzz**2
  NrotN = - Nr*Nfz + Nf*(Nrz-Nzr) + Nz*(Nf/r + Nfr)
  Fg = XiH**2 * (gN2 - 1D0/16D0*(sqrt(3D0)*DivN + sqrt(5D0)*NRotN )**2)

  !! Full energy
  F = Fh + Fg
!  F=gN2
!  F  =  (Nrr)**2 + (Nfr)**2 + Nzr**2

end

!!! User defined texture limits and boundary conditions
subroutine VAR_LIM(r,z,n, Umin,Umax)
  implicit none
  integer n
  double precision r,z,Umin,Umax

  double precision pi
  data pi/3.14159265358979/

  !! limits and boundary conditions
  !! for \alpha (n==1) and \beta (n==2)

  if (n.eq.1) then
    Umin = -PI
    Umax = PI

    if (r.ge.1D0) then
      Umin = acos(0.5D0)
      Umax = Umin
    endif

!    if (z.ge.1D0) then
!      Umin = 0D0
!      Umax = Umin
!    endif

  else
    Umin = 0
    Umax = PI
!    if (z.ge.1D0) then
!      Umin = 0D0
!      Umax = Umin
!    endif
    if (r.le.1D-7) then
      Umin = 0D0
      Umax = Umin
    endif
    if (r.ge.1D0) then
      Umin = acos(1D0/sqrt(5D0))
      Umax = Umin
    endif
  endif
end