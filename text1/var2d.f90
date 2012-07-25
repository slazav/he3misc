!!! Solve variation problem with fixed step H, and number of points
!  U     -- on input: initial values, on output: result
!  NU    -- U-function dimensions

subroutine VAR2D(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,HS,HF)
  implicit none
  integer NX,NY,NU, NXX, NYY, N, X,Y
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX,NY,NU), HS,HF,H

!  NXX=1
!  NYY=1
!  do while (NXX.lt.NX.and.NYY.lt.NY)
!    NXX=NXX*2
!    NYY=NYY*2
!    write(*,*) 'Nx,Ny', NXX, NYY
!    if (NXX.gt.NX) NXX=NX
!    if (NYY.gt.NY) NYY=NY
    ! interpolate U to new grid
    !...
    ! set new HS
    call VAR2D_PASS2(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,HS,HF)
!  enddo
end

!!! Do minimization at fixed grid
subroutine VAR2D_PASS2(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,HS,HF)
  implicit none
  integer NX,NY,NU
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX,NY,NU), HS,HF,H

  H=HS*2D0
  do while (H.gt.HF)
    H=H/4D0
    if (H.lt.HF) H=HF
    call VAR2D_PASS1(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,H)
  enddo
end

!!! Do single pass through all data at fixed NX,NY and H
subroutine VAR2D_PASS1(XMIN,XMAX,YMIN,YMAX,U,NX,NY,NU,H)
  implicit none
  integer NX,NY,NU, IX,IY,IU, N
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX,NY,NU), H, DX,DY

  DX = (XMAX-XMIN)/(NX-1) !! we need 1/2 DX space on both sides of Xmin..Xmax
  DY = (YMAX-YMIN)/(NY-1)

  N=1
  do while (N.gt.0)
    N=0
    do IU = 1,NU
      do IX = 1,NX
        do IY = 1,NY
          N = N + SETU(IX, IY, IU, H)
        enddo
      enddo
    enddo
  enddo


  contains


  !! calculate and compare integrals in 2x2 cell area for
  !! -H,0,+H steps and set new value for U
  !! returns 1 if value was changed
  function SETU(IX, IY, IU, H)
    implicit none
    integer IX,IY,IU,SETU
    double precision H

    double precision S0, SP, SM
    double precision U_, Umin, Umax
    double precision X_, Y_
    X_ = XMIN + DX * (IX-1)
    Y_ = YMIN + DY * (IY-1)
    if (IX.eq.NX) X_=XMAX  ! exact values for boundary conditions
    if (IY.eq.NY) Y_=YMAX

    call VAR_LIM(X_,Y_,IU, Umin,Umax) ! get min/max values

    ! fix value if needed
    if (U(IX, IY, IU).lt.Umin) U(IX, IY, IU) = Umin
    if (U(IX, IY, IU).gt.Umax) U(IX, IY, IU) = Umax

    U_ = U(IX, IY, IU) ! save initial value

    S0 = INT4(IX, IY)

    SETU=1 ! return 1 by default

    if (U_ + H.le.Umax) then
      U(IX, IY, IU) = U_ + H
      SP = INT4(IX,IY)
      if (SP.lt.S0) return  ! U is already set to U_ + H
    endif

    if (U_ - H.ge.Umin) then
      U(IX, IY, IU) = U_ - H
      SM = INT4(IX,IY)
      if (SM.lt.S0) return  ! U is already set to U_ - H
    endif

    ! set unmodified value, return 0
    U(IX, IY, IU) = U_
    SETU=0
    return
  end

  !! calculate integral in 2x2 cells
  !! which includes IX,IY function value
  function INT4(IX, IY)
    implicit none
    integer IX,IY, IIX, IIY
    double precision INT4

    INT4=0
    do IIX = max(1,IX-1),min(IX,NX-1)
      do IIY = max(1,IY-1),min(IY,NY-1)
        INT4 = INT4 + INT1(IIX,IIY)
      enddo
    enddo
  end


  !! calculate integral in one cell IX..IX+1, IY..IY+1
  function INT1(IX, IY)
    implicit none
    integer IX,IY, IU
    double precision INT1

    double precision U_(NU), UX_(NU), UY_(NU)
    double precision X_, Y_

    ! values in the center of the cell:
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
