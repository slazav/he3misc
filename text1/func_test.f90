!!! User-defined functional for minimization
!
! 2D-texture energies in rotating 3He-B
! as a function of alpha and beta angles of n-vector
! (Gradient energy +  Orientation energy)
subroutine VAR_FNC(x,y,U,Ux,Uy, NU, F)
  implicit none
  integer NU
  double precision x,y,U(NU),Ux(NU),Uy(NU),F

  F = Ux(1)**2 + Uy(1)**2 + Ux(2)**2 + Uy(2)**2
end

!!! User defined texture limits and boundary conditions
subroutine VAR_LIM(x,y,n, Umin,Umax)
  implicit none
  integer n
  double precision x,y,Umin,Umax

  Umin = -1
  Umax = 1

  if (x.ge.1D0.or.x.le.0D0.or.y.ge.1D0.or.y.le.0D0) then
    Umin = 0.1D0
    Umax = Umin
  endif

  Umin = 1 - 16D0*(x-0.5D0)**2 - 16D0*(y-0.5D0)**2
!  Umax = 16D0*(x-0.5D0)**2 + 16D0*(y-0.5D0)**2 - 1

end
