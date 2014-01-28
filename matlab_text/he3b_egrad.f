!     Compare he3b_egrad_nt3d and he3b_egrad_n2d
      program test1
        real*8 X(4), DX(4,3), L1, L2
        real*8 E1, EX1(4), EDX1(4,3)
        real*8 E2, EX2(4), EDX2(4,3)
        real*8 tol
        integer i1,i2,i3,i4,i5,i6,i7,i8,i9
        L1=1.1
        L2=3.2

!       Theta=const
        X(4) = acos(-0.25D0)
        DX(4,1) = 0D0
        DX(4,2) = 0D0
        DX(4,3) = 0D0
!       d/dz=0
        DX(1,3) = 0D0
        DX(2,3) = 0D0
        DX(3,3) = 0D0

        call srand(0)
        tol=1D-10

        do i = 1,1 00
          X(1) = rand(0)*2D0-1D0
          X(2) = rand(0)*2D0-1D0
          X(3) = rand(0)*2D0-1D0
          DX(1,1) = rand(0)*2D0-1D0
          DX(2,1) = rand(0)*2D0-1D0
          DX(3,1) = rand(0)*2D0-1D0
          DX(1,2) = rand(0)*2D0-1D0
          DX(2,2) = rand(0)*2D0-1D0
          DX(3,2) = rand(0)*2D0-1D0
          call he3b_egrad_n2d (L1, L2, X, DX, E1, EX1, EDX1)
          call he3b_egrad_nt3d(L1, L2, X, DX, E2, EX2, EDX2)
          if (dabs(E1-E2).gt.tol)
     .       write(*,*) 'E>> ', X(1), X(2), X(3), E1-E2
          if (dabs(EX1(1)-EX2(1)).gt.tol)
     .       write(*,*) 'EX1>> ', X(1), X(2), X(3), EX1(1)-EX2(1)
          if (dabs(EX1(2)-EX2(2)).gt.tol)
     .       write(*,*) 'EX1>> ', X(1), X(2), X(3), EX1(2)-EX2(2)
          if (dabs(EX1(3)-EX2(3)).gt.tol)
     .       write(*,*) 'EX1>> ', X(1), X(2), X(3), EX1(3)-EX2(3)
        enddo
      end

! He3-B gradient energy as a function of nx,ny,nz,theta.
! in 3D coordinates x,y,z
!
! E_G = L1 dR_ki/dx_i dR_kj/dx_j + L2 dR_kj/dx_i dR_ki/dx_j
! where R is a rotation matrix R(nx, ny, nz, theta).
! Input:
!   X(1:4)      - nx,ny,nz,theta
!   DX(1:4,1:3) - dX_i/dx_j = d(nx,ny,nz,theta)/d(x,y,z)
!   L1, L2      - coefficients
! Output:
!   E            - energy
!   EX(1:4)      - derivatives dE/dXi
!   EDX(1:4,1:3) - derivatives dE/d(DXij)
! There is a PDF with formulas somewhere
! V.Zavjalov 01.2014
! calculation time: 1.0e-4 s on rota computer

      subroutine he3b_egrad_nt3d(L1, L2, X, DX, E, EX, EDX)
        real*8 X(4), DX(4,3), L1, L2
        real*8 E, EX(4), EDX(4,3)

        integer i,j,k, a,b,c ! loop variables
        real*8 ct,st,v1,v2,v3

!       dd is Kronecker delta, ee is Levi-Civita symbol,
!       en = e_ijk n_k
!       DR=dR/dX, D2R=d2R/dX2 - no need to initialize.
        real*8 ee(3,3,3), dd(3,3), en(3,3), DR(3,3,4), D2R(3,3,4,4)
        data ee(1,2,3), ee(2,3,1), ee(3,1,2) /3*1D0/,
     .       ee(3,2,1), ee(2,1,3), ee(1,3,2) /3*-1D0/,
     .      (dd(i,i),i=1,3) /3*1D0/,
     .    ee /27*0D0/, dd /9*0D0/

        ct=cos(X(4));
        st=sin(X(4));

!       calculate en matrix (3x3 = 9 loops)
        do i=1,3
          do j=1,3
            en(i,j) = 0
            do a=1,4
              en(i,j) = en(i,j) + ee(i,j,a)*X(a)
            enddo
          enddo
        enddo

!       calculate DR matrix (3x3x4 = 32 loops)
        do i=1,3
          do j=1,3
            do a=1,4
              if (a.ne.4) then
                DR(i,j,a) = (1D0-ct)*(dd(i,a)*X(j) + dd(j,a)*X(i))
     .             - st*ee(i,j,a)
              else
                DR(i,j,a) = st*(X(i)*X(j)-dd(i,j)) - ct*en(i,j)
              endif
            enddo
          enddo
        enddo

!       calculate D2R matrix (3x3x4x4 = 144 loops)
        do i=1,3
          do j=1,3
            do a=1,4
              do b=1,4
                if (a.ne.4.and.b.ne.4) then
                  D2R(i,j,a,b) =
     .              (1D0-ct)*(dd(i,a)*dd(j,b)+dd(j,a)*dd(i,b))
                elseif (a.ne.4.and.b.eq.4) then
                  D2R(i,j,a,b) =
     .               st*(dd(i,a)*X(j)+dd(j,a)*X(i)) - ct*ee(i,j,a)
                elseif (a.eq.4.and.b.ne.4) then
                  D2R(i,j,a,b) =
     .               st*(dd(i,b)*X(j)+dd(j,b)*X(i)) - ct*ee(i,j,b)
                else
                  D2R(i,j,a,b) =
     .               ct*(X(i)*X(j)-dd(i,j)) - ct*en(i,j)
                endif
              enddo
            enddo
          enddo
        enddo

!	initialize EX - it's better to do it here
        do a=1,4
          EX(a)=0D0
        enddo

!       calculate result (4x3x4x3x3x4 = 1728 loops)
        E=0D0
        do a=1,4
          do i=1,3
            EDX(a,i)=0D0
            do b=1,4
              do j=1,3
                v1=DX(a,i)*DX(b,j)
                v2=2D0*DX(b,j)
                do k=1,3
                  v3 = L1*DR(k,i,a)*DR(k,j,b) + L2*DR(k,j,a)*DR(k,i,b)
                  E = E + v1*v3
                  EDX(a,i) = EDX(a,i) + v2*v3
                  do c=1,4
                    EX(c) = v1*(
     .                L1*(D2R(k,i,a,c)*DR(k,j,b)
     .                  + D2R(k,j,b,c)*DR(k,i,a))
     .              + L2*(D2R(k,j,b,c)*DR(k,i,b)
     .                  + D2R(k,j,a,c)*DR(k,j,a)))
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

      end

! He3-B gradient energy as a function of nx,ny,nz (theta=acos(-1/4)).
! in 2D coordinates x,y.
!
! E_G = L1 dR_ki/dx_i dR_kj/dx_j + L2 dR_kj/dx_i dR_ki/dx_j
! where R is a rotation matrix R(nx, ny, nz, theta).
! Input:
!   n(1:3)      - nx,ny,nz
!   dn(1:3,1:2) - d(nx,ny,nz)/d(x,y)
!   L1, L2      - coefficients
! Output:
!   e            - energy
!   en(1:3)      - derivatives de/dn_i
!   edn(1:3,1:2) - derivatives de/d(dn_ij)
! There is a PDF with formulas somewhere
! V.Zavjalov 01.2014
! calculation time: 3.32e-6 s on rota computer

      subroutine he3b_egrad_n2d(l1, l2, n,dn, e, en,edn)
        real*8 l1,l2,n(3),dn(3,2), e, en(3), edn(3,2)
        real*8 c1,c2,c3,c4,c5
        c1 = 25D0/16D0*(l1+l2)
        c2 = -5D0*sqrt(5D0)/4D0*(l1+l2)
        c3 = -5D0*sqrt(5D0)/4D0
        c4 = 5D0/16D0*(l1+l2)
        c5 = 5D0/4D0

        e = c1 * ( (n(1)*dn(1,1) + n(2)*dn(1,2))**2                      &
     &           + (n(1)*dn(2,1) + n(2)*dn(2,2))**2                      &
     &           + (n(1)*dn(3,1) + n(2)*dn(3,2))**2 )                    &
     &    + c2 * (n(1)*dn(3,2)*dn(2,2) - n(2)*dn(3,1)*dn(1,1)            &
     &          + n(3)*(dn(1,1)+dn(2,2))*(dn(2,1)-dn(1,2)))              &
     &    + c3 * ( l1*(n(1)*dn(3,2)*dn(1,1) - n(2)*dn(3,1)*dn(2,2))      &
     &           + l2*(n(1)*dn(3,1)*dn(1,2) - n(2)*dn(3,2)*dn(2,1)) )    &
     &    + c4 * ( 5D0*dn(1,1)*dn(1,1) + 5D0*dn(2,2)*dn(2,2)             &
     &           + 3D0*dn(2,1)*dn(2,1) + 3D0*dn(1,2)*dn(1,2) )           &
     &    + c5 * ( (5D0*l1-3D0*l2) * dn(1,1)*dn(2,2)                     &
     &           + (5D0*l2-3D0*l2) * dn(1,2)*dn(2,1) )

        en(1) = 2D0*c1 *                                                 &
     &     ( (n(1)*dn(1,1) + n(2)*dn(1,2))*dn(1,1)                       &
     &     + (n(1)*dn(2,1) + n(2)*dn(2,2))*dn(2,1)                       &
     &     + (n(1)*dn(3,1) + n(2)*dn(3,2))*dn(3,1) )                     &
     &    + c2 * dn(3,2)*dn(2,2)                                         &
     &    + c3 * ( l1*dn(3,2)*dn(1,1) + l2*dn(3,1)*dn(1,2) )
        en(2) = 2D0*c1 *                                                 &
     &     ( (n(1)*dn(1,1) + n(2)*dn(1,2))*dn(1,2)                       &
     &     + (n(1)*dn(2,1) + n(2)*dn(2,2))*dn(2,2)                       &
     &     + (n(1)*dn(3,1) + n(2)*dn(3,2))*dn(3,2) )                     &
     &    - c2 * dn(3,1)*dn(1,1)                                         &
     &    - c3 * ( l1*dn(3,1)*dn(2,2) + l2*dn(3,2)*dn(2,1) )
        en(3) = c2 * (dn(1,1)+dn(2,2))*(dn(2,1)-dn(1,2))

        edn(1,1) = 2D0*c1 * (n(1)*dn(1,1) + n(2)*dn(1,2))*n(1)           &
     &    + c2 * (- n(2)*dn(3,1) + n(3)*(dn(2,1)-dn(1,2)))               &
     &    + c3 * l1*n(1)*dn(3,2)                                         &
     &    + c4 * 10D0*dn(1,1)*dn(1,1)                                    &
     &    + c5 * (5D0*l1-3D0*l2) * dn(2,2)
        edn(1,2) = 2D0*c1 * (n(1)*dn(1,1) + n(2)*dn(1,2))*n(2)           &
     &    - c2 * n(3)*(dn(1,1)+dn(2,2))                                  &
     &    + c3 * l2*n(1)*dn(3,1)                                         &
     &    + c4 * 6D0*dn(1,2)                                             &
     &    + c5 * (5D0*l2-3D0*l2) * dn(2,1)
        edn(2,1) = 2D0*c1 * (n(1)*dn(2,1) + n(2)*dn(2,2))*n(1)           &
     &    + c2 * n(3)*(dn(1,1)+dn(2,2))                                  &
     &    - c3 * l2*n(2)*dn(3,2)                                         &
     &    + c4 * 6D0*dn(2,1)                                             &
     &    + c5 * (5D0*l2-3D0*l2) * dn(1,2)
        edn(2,2) = 2D0*c1 * (n(1)*dn(2,1) + n(2)*dn(2,2))*n(2)           &
     &    + c2 * (n(1)*dn(3,2) + n(3)*(dn(2,1)-dn(1,2)))                 &
     &    - c3 * n(2)*dn(3,1)                                            &
     &    + c4 * 10D0*dn(2,2)                                            &
     &    + c5 * (5D0*l1-3D0*l2) * dn(1,1)
        edn(3,1) = 2D0*c1 * (n(1)*dn(3,1) + n(2)*dn(3,2))*n(1)           &
     &    - c2 * n(2)*dn(1,1)                                            &
     &    - c3 * ( l1*n(2)*dn(2,2) - l2*n(1)*dn(1,2))
        edn(3,2) = 2D0*c1 * (n(1)*dn(3,1) + n(2)*dn(3,2))*n(2)           &
     &    + c2 * n(1)*dn(2,2)                                            &
     &    + c3 * ( l1*n(1)*dn(1,1) - l2*n(2)*dn(2,1))
      end

