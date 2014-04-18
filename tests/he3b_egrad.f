      program test1
        real*8 X(4), DX(4,3), L1, L2
        real*8 E, EX(4), EDX(4,3)
        integer i
        do i=1,1000000
          call he3b_egrad(X, DX, L1, L2,  E, EX, EDX)
        enddo
      end

! He3-B gradient energy
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

      subroutine he3b_egrad(X, DX, L1, L2,  E, EX, EDX)
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
