!!! Solve variation problem with fixed step H, and number of points
!  U     -- on input: initial values, on output: result
!  NU    -- U-function dimensions
! User must provide two subroutines:
! VAR_FNC(X,Y,U,UY,UY, NU, F)
! VAR_LIM(X,Y,IU, UMIN,UMAX) -- calculate limits and soundary conds UMIN,UMAX for X,Y,IU
!

subroutine VAR2D(XMIN,XMAX,YMIN,YMAX,U,NX0,NY0,NU,HF)
  implicit none
  integer NX0,NY0,NX,NY,NU, NXP,NYP, N, X,Y
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX0,NY0,NU), HS,HF,H, UDIFF

  NX=3
  NY=3
  do while (NX.lt.NX0.and.NY.lt.NY0)
    NXP=NX ! save previous values
    NYP=NY
    NX=2*NX-1 ! new grid
    NY=2*NY-1
    if (NX.gt.NX0) NX=NX0
    if (NY.gt.NY0) NY=NY0

    HS = UDIFF(U, NX0,NY0,NXP,NYP,NU)/2D0
    call INTERP(U, NX0,NY0,NU, NXP,NYP, NX,NY)

    if (HS.le.0) HS=1D0
    write(*,*) ' NX,NY= ', NX, NY
    call VAR2D_PASS2(XMIN,XMAX,YMIN,YMAX,U,NX0,NY0,NX,NY,NU,HS,HF)
  enddo

end

! get max function difference
function UDIFF(U, NX0,NY0,NX,NY,NU)
  implicit none
  double precision U(NX0,NY0,NU), UDIFF
  integer NX0,NY0,NX,NY,NU, IU,IX,IY
  double precision DD
  UDIFF=0
  do IU = 1,NU
    do IX = 2,NX-1
      do IY = 2,NY-1
        DD = max( &
          abs(2*U(IX,IY,IU) - U(IX+1,IY,IU) - U(IX-1,IY,IU)), &
          abs(2*U(IX,IY,IU) - U(IX,IY+1,IU) - U(IX,IY-1,IU)))
        if (UDIFF.lt.DD) UDIFF=DD
      enddo
    enddo
  enddo
end

! interpolate U(NX,NY,NU) from NX1 x NY1 to NX2 x NX2.
! NX1 <= NX2 <= NX,  NY1 <= NY2 <= NY 
subroutine INTERP(U, NX,NY,NU, NX1,NY1, NX2,NY2)
  implicit none
  double precision U(NX,NY,NU)
  integer NX,NY,NU, NX1,NY1, NX2,NY2

  integer IX,IY,IU, IX1,IY1
  double precision SX,SY ! scale
  double precision DX,DY, U00, U10, U01, U11, U1, U2

  SX = (NX1-1D0)/(NX2-1D0)
  SY = (NY1-1D0)/(NY2-1D0)

  do IU = 1,NU
    do IX = NX2,1,-1
      do IY = NY2,1,-1
        IX1 = int(floor((IX-1)*SX))+1
        IY1 = int(floor((IY-1)*SY))+1
        DX = (IX-1)*SX+1D0-IX1
        DY = (IY-1)*SY+1D0-IY1

        U00 = U(IX1,IY1,IU)
        U10 = U(IX1+1,IY1,IU)
        U01 = U(IX1,IY1+1,IU)
        U11 = U(IX1+1,IY1+1,IU)
        U1 = U00 + DX*(U10-U00)
        U2 = U01 + DX*(U11-U01)

        U(IX,IY,IU) = U1 + (U2-U1)*DY
      enddo
    enddo
  enddo
end

!!! Do minimization at fixed grid
subroutine VAR2D_PASS2(XMIN,XMAX,YMIN,YMAX,U,NX0,NY0,NX,NY,NU,HS,HF)
  implicit none
  integer NX0,NY0,NX,NY,NU
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX0,NY0,NU), HS,HF,H

  H=HS*2D0
  do while (H.gt.HF)
    H=H/4D0
    if (H.lt.HF) H=HF
    write(*,*) ' H= ', H
    call VAR2D_PASS1(XMIN,XMAX,YMIN,YMAX,U,NX0,NY0,NX,NY,NU,H)
  enddo
end

!!! Do single pass through all data at fixed NX,NY and H
!!! NX0, NY0, NU -- U dimensions
!!! NX, NY -- current grid size
subroutine VAR2D_PASS1(XMIN,XMAX,YMIN,YMAX,U,NX0,NY0,NX,NY,NU,H)
  implicit none
  integer NX0,NY0,NX,NY,NU, IX,IY,IU, N
  double precision XMIN,XMAX,YMIN,YMAX
  double precision U(NX0,NY0,NU), H, H1, DX,DY

  DX = (XMAX-XMIN)/(NX-1)
  DY = (YMAX-YMIN)/(NY-1)

  N=1
  do while (N.gt.0)
    N=0
    H1=H
    do IU = 1,NU
      do IX = 1,NX
        do IY = 1,NY
          H1 = SETU(IX, IY, IU, H1)
          if (H1.ne.0D0) then
            N = N + 1
          else
            H1=H
          endif
        enddo
      enddo
    enddo
  enddo

  contains

  !! calculate and compare integrals in 2x2 cell area for
  !! -H,0,+H steps and set new value for U
  !! returns step used
  function SETU(IX, IY, IU, H)
    implicit none
    integer IX,IY,IU
    double precision H, SETU

    double precision S0, SP, SM
    double precision U_, Umin, Umax
    double precision X_, Y_

    X_ = XMIN + DX * (IX-1)
    Y_ = YMIN + DY * (IY-1)
    if (IX.eq.NX) X_=XMAX  ! exact values for boundary conditions
    if (IY.eq.NY) Y_=YMAX

    call VAR_LIM(X_,Y_,IU, Umin,Umax) ! get min/max values

    ! fix value if needed
    SETU = 0D0
    if (U(IX, IY, IU).lt.Umin) then
      U(IX, IY, IU) = Umin
      return
    endif
    if (U(IX, IY, IU).gt.Umax) then
      U(IX, IY, IU) = Umax
      return
    endif

    U_ = U(IX, IY, IU) ! save initial value
    S0 = INT4(IX,IY)

    if (U_+H.le.Umax.and.U_+H.ge.Umin) then
      U(IX, IY, IU) = U_ + H
      SP = INT4(IX,IY)
      if (SP.lt.S0) then
         SETU = H
         return  ! U is already set to U_ + H
      endif
    endif

    if (U_-H.le.Umax.and.U_-H.ge.Umin) then
      U(IX, IY, IU) = U_ - H
      SM = INT4(IX,IY)
      if (SM.lt.S0) then
        SETU = -H
        return  ! U is already set to U_ - H
      endif
    endif

    ! set unmodified value, return 0
    U(IX, IY, IU) = U_
    SETU=0D0
    return
  end

  !! calculate integral in 2x2 cells
  !! which includes IX,IY function value
  function INT4(IX, IY)
    implicit none
    integer IX,IY, IIX, IIY
    double precision INT4

    INT4=0
    if (IX.ne.NX) then
      if (IY.ne.NY) then
        INT4 = INT4 + INT1(IX,IY,    1, 1)
        INT4 = INT4 + INT1(IX+1,IY, -1, 1)
        INT4 = INT4 + INT1(IX,IY+1,  1,-1)
      endif
      if (IY.ne.1) then
        INT4 = INT4 + INT1(IX,IY,    1,-1)
        INT4 = INT4 + INT1(IX+1,IY, -1,-1)
        INT4 = INT4 + INT1(IX,IY-1,  1, 1)
      endif
    endif
    if (IX.ne.1) then
      if (IY.ne.NY) then
        INT4 = INT4 + INT1(IX,IY,   -1, 1)
        INT4 = INT4 + INT1(IX-1,IY,  1, 1)
        INT4 = INT4 + INT1(IX,IY+1, -1,-1)
      endif
      if (IY.ne.1) then
        INT4 = INT4 + INT1(IX,IY,   -1,-1)
        INT4 = INT4 + INT1(IX-1,IY,  1,-1)
        INT4 = INT4 + INT1(IX,IY-1, -1, 1)
      endif
    endif
  end


  !! calculate integral in one cell IX..IX+DIX, IY..IY+DIY
  function INT1(IX, IY, DIX, DIY)
    implicit none
    integer IX,IY, DIX,DIY, I
    double precision INT1

    double precision U_(NU), UX_(NU), UY_(NU)
    double precision X_, Y_

    ! values in the center of the cell:
    X_ = XMIN + DX * (IX-1D0)
    Y_ = YMIN + DY * (IY-1D0)
    do I=1,NU
      U_(I)  = U(IX,IY,I)
      UX_(I) = (U(IX+DIX,IY,I) - U(IX,IY,I))/DIX/DX
      UY_(I) = (U(IX,IY+DIY,I) - U(IX,IY,I))/DIY/DY
    enddo
    call VAR_FNC(X_, Y_, U_, UX_, UY_, NU, INT1)
!    INT1 = DX*DY*INT1
  end

end
