program texture
  implicit none

  common /texture_pars/ XiH, Lambda, KappaH
  double precision XiH, Lambda, KappaH
  data XiH    /1D0/    &
       Lambda /1D0/ &
       KappaH /1D0/

  integer nr,nz, ir,iz
  parameter (nr=100, nz=100)

  double precision U(nr,nz,2), H
  double precision rmin,rmax,zmin,zmax
  data rmin /0D0/ &
       rmax /1D0/ &
       zmin /0D0/ &
       zmax /1D0/ &
       H/0.01/

  do ir=1,nr
    do iz=1,nz
      U(ir,iz,1) = 0D0;
      U(ir,iz,2) = 0D0;
    enddo
  enddo

  call VAR2D(rmin,rmax,zmin,zmax,U,nr,nz,2,H)

  do ir=1,nr
    do iz=1,nz
      write(*,*) ir, iz, U(ir,iz,1), U(ir,iz,2)
    enddo
  enddo

end

!!! Solve variation problem with fixed step H, and number of points
!  U     -- on input: initial values, on output: result
!  NU    -- U-function dimensions
subroutine VAR2D(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,H)
  implicit none
  integer NX,NY,NU, IX,IY,IU, N
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX,NY,NU), H, DX,DY

  DX = (XMAX-XMIN)/(NX-1)
  DY = (YMAX-YMIN)/(NY-1)

  N=1
  do while (N.gt.0)
    N=0
    do IU = 1,NU
      do IX = 1,NX-1
        do IY = 1,NY-1
          N = N + SETU(IX, IY, IU)
        enddo
      enddo
    enddo
  enddo


  contains


  !! calculate and compare integrals in 3x3 cell area for
  !! -H,0,+H steps and set new value for U
  !! returns 1 if value was changed
  function SETU(IX, IY, IU)
    implicit none
    integer IX,IY,IU,SETU

    double precision S0, SP, SM
    double precision U_, Umin, Umax
    double precision X_, Y_

    X_ = XMIN + DX * (IX-0.5D0)
    Y_ = YMIN + DY * (IY-0.5D0)

    call VAR_LIM(X_,Y_,IU, Umin,Umax) ! get min/max values

    ! fix value
    if (U(IX, IY, IU).lt.Umin) U(IX, IY, IU) = Umin
    if (U(IX, IY, IU).gt.Umax) U(IX, IY, IU) = Umax

    U_ = U(IX, IY, IU) ! save initial value

    S0 = INT9(IX, IY)
    SETU=1

    if (U_ + H.le.Umax) then
      U(IX, IY, IU) = U_ + H
      SP = INT9(IX,IY)
    else
      SP = S0+1D0 ! something greater than S0
    endif

    if (SP.lt.S0) return  ! U is already set to U_ + H

    if (U_ - H.ge.Umin) then
      U(IX, IY, IU) = U_ - H
      SM = INT9(IX,IY)
    else
      SM = S0+1D0
    endif

    if (SM.lt.S0) return  ! U is already set to U_ - H

    U(IX, IY, IU) = U_
    SETU=0
    return
  end

  !! calculate integral in 3x3 cell area
  function INT9(IX, IY)
    implicit none
    integer IX,IY, IIX, IIY
    double precision INT9

    INT9=0
    do IIX = max(1,IX-1),min(NX,IX+1)
      do IIY = max(1,IY-1),min(NY,IY+1)
        INT9 = INT9 + INT1(IIX,IIY)
      enddo
    enddo
  end


  !! calculate integral in one cell
  function INT1(IX, IY)
    implicit none
    integer IX,IY, IU
    double precision INT1

    double precision U_(NU), UX_(NU), UY_(NU)
    double precision X_, Y_

    ! values in the center of the cell IX..IX+1, IY..IY+1:
    X_ = XMIN + DX * (IX-0.5D0)
    Y_ = YMIN + DY * (IY-0.5D0)
    do IU=1,NU
      U_(IU)  = ( U(IX,IY,IU) + U(IX+1,IY,IU) + &
                  U(IX,IY+1,IU) + U(IX+1,IY+1,IU))/4D0
      UX_(IU) = ( U(IX+1,IY+1,IU) + U(IX+1,IY,IU) &
                - U(IX,IY+1,IU) - U(IX,IY,IU))/2D0/DX
      UY_(IU) = ( U(IX+1,IY+1,IU) + U(IX,IY+1,IU) &
                - U(IX+1,IY,IU) - U(IX,IY,IU))/2D0/DY
    enddo

    call VAR_FNC(X_, Y_, U_, UX_, UY_, NU, INT1)
    INT1 = DX*DY*INT1
  end

end


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
  double precision Nr,Nf,Nz, Nrr,Nfr,Nzr, Nrz,Nfz,Nzz
  double precision NxZ2, divN, NrotN, Fh, Fg

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
    !! boundary conditions:
    if (r.ge.1D0) then
      Umin = PI/3D0
      Umax = Umin
    endif
  else
    Umin = 0D0
    Umax = 2D0*PI
    !! boundary conditions
    if (r.le.0D0) then
      Umin = 0D0
      Umax = Umin
    endif
    if (r.ge.1D0) then
      Umin = acos(1D0/sqrt(5D0))
      Umax = Umin
    endif
  endif
end

