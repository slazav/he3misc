C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
C       T - time
C       X - x-coord
C       U   - Mx My Mz Nx Ny Nz Theta
C       UX  - dU/dx
C       UXX - d2U/dx2
C       FV  - result (dU/dt)
        include 'vmcw.fh'
        include 'he3.fh'
        include 'legg_eq/he3b_legg_rot1d.fh'
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)

C       set parameters
        Hx = HR0+HR_SWR*T
        Hz = H + GRAD*X
        W0 = he3_gyro*(H + GRAD*(LP0+LP_SWR*T))
        Flegg = LF0 + LF_SWR*T
        Diff  = DF0 + DF_SWR*T
        Cpar  = CPAR0 + CPAR_SWR*T
        dCpar = 0
        Tf    = TF0+TF_SWR*T
        T1    = T11

C       spatial modulation
        Cpar0=Cpar
        dCpar = -0.25D0 * Cpar0 * AER_STEP(X,1)/
     .    dsqrt(1D0 - 0.5D0 * AER_STEP(X,0))
        Cpar  = Cpar0 * dsqrt(1D0 - 0.5D0 * AER_STEP(X,0))
        Flegg = Flegg * dsqrt(1D0 - 0.835D0 * AER_STEP(X,0))
        Tf    = Tf * (1D0 - 0.5D0 * AER_STEP(X,0))

        call he3b_legg_rot1d(U,UX,UXX,FV)

      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        include 'vmcw.fh'
        include 'he3.fh'

        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

        if(IBN.EQ.2)THEN       ! CLOSED CELL

C         fix n vector length
          UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UNx=U(4)/UN
          UNy=U(5)/UN
          UNz=U(6)/UN

          ST=dsin(U(7))
          ST2=2.0D0*ST
          CT=dcos(U(7))
          CTM=1.0D0-CT
          CTM2=2.0D0*CTM
          DD45=UNx*UX(5)-UX(4)*UNy
          FTN=CTM*DD45-ST*UX(6)-UX(7)*UNz
          CTF=CTM*FTN
          STF=ST*FTN
          FTN4=CTM*UX(5)
          FTN5=-CTM*UX(4)
          FTN7=ST*DD45-CT*UX(6)
          FTNX4=-CTM*UNy
          FTNX5=CTM*UNx
          C46=CTM*UNx*UNz+UNy*ST
          C56=CTM*UNy*UNz-UNx*ST             !!!!!!!!!!!
          C66=CTM*UNz**2+CT
          C266=2.0D0-C66

          W=He3_gyro*H
!          W=He3_gyro*(H + GRAD*X)
          AF=-(CPAR0+CPAR_SWR*T)**2/W
          AF=AF*(1D0 - 0.5D0 * AER_STEP(X,0))
          DA=-(DF0+DF_SWR*T)/AF

          DBDUX(4,1)=DA
          DBDUX(5,2)=DA
          DBDUX(6,3)=DA

          DBDU(4,4)=2.0D0*UX(7)+CTF*UNz+C46*FTN4
          DBDU(4,5)=CTM2*UX(6)+STF+C46*FTN5
          DBDU(4,6)=-CTM2*UX(5)+CTF*UNx-C46*UX(7)
          DBDU(4,7)=2.0D0*(CT*UX(4)+ST*(UNy*UX(6)-UX(5)*UNz))+
     *     STF*UNx*UNz+UNy*CT*FTN+C46*FTN7

          DBDU(5,4)=-CTM2*UX(6)-STF+C56*FTN4
          DBDU(5,5)=2.0D0*UX(7)+CTF*UNz+C56*FTN5
          DBDU(5,6)=CTM2*UX(4)+CTF*UNy-C56*UX(7)
          DBDU(5,7)=2.0D0*(CT*UX(5)-ST*(UNx*UX(6)-UX(4)*UNz))+
     *     STF*UNy*UNz-UNx*CT*FTN+C56*FTN7

          DBDU(6,4)=CTM2*UX(5)+C66*FTN4
          DBDU(6,5)=-CTM2*UX(4)+C66*FTN5
          DBDU(6,6)=2.0D0*UNz*CTF+C266*UX(7)
          DBDU(6,7)=2.0D0*(CT*UX(6)+ST*DD45)+STF*(UNz**2-1.0D0)+C66*FTN7

          DBDUX(4,4)=ST2+C46*FTNX4
          DBDUX(4,5)=-CTM2*UNz+C46*FTNX5
          DBDUX(4,6)=CTM2*UNy-C46*ST
          DBDUX(4,7)=2.0D0*UNx-C46*UNz

          DBDUX(5,4)=CTM2*UNz+C56*FTNX4
          DBDUX(5,5)=ST2+C56*FTNX5
          DBDUX(5,6)=-CTM2*UNx-C56*ST
          DBDUX(5,7)=2.0D0*UNy-C56*UNz

          DBDUX(6,4)=-CTM2*UNy+C66*FTNX4
          DBDUX(6,5)=CTM2*UNx+C66*FTNX5
          DBDUX(6,6)=C266*ST
          DBDUX(6,7)=C266*UNz
C          DBDU(7,4)=UX(4)         !!
C          DBDU(7,5)=UX(5)         !!
C          DBDU(7,6)=UX(6)         !!
C          DBDUX(7,4)=UNx         !!
C          DBDUX(7,5)=UNy         !!
C          DBDUX(7,6)=UNz         !!
        endif
        return
      end
C-- SET_ICOND -- INITIAL CONDITIONS ---------------------------------
      subroutine SET_ICOND()
        include 'vmcw.fh'
        include 'par.fh'
        include 'he3.fh'
        real*8 USOL, XSOL
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)

        real*8 BET, DELTA, DELTAX, DELTAY,
     .         UCTG, UNX,UNY,UNZ, UMX,UMY,UMZ,
     .         UNZ2
        integer I,J,K

        BET=BETA*PI/180.0D0
        UMZ=dcos(BET)
        UMX=dsin(BET)*dsqrt(0.5D0)
        UMY=UMX
        if(UMZ.GE.-0.25D0)THEN
          UNZ2=0.8D0*(0.25D0+dcos(BET))
          UNZ=dsqrt(UNZ2)
          DELTA=(25.0D0*UNZ2+15.0D0)/16.0D0
          DELTAX=dsin(BET)*dsqrt(0.5D0)*
     *     (UNZ*1.25D0-dsqrt(15.0D0)*0.25D0)
          DELTAY=dsin(BET)*dsqrt(0.5D0)*
     *     (UNZ*1.25D0+dsqrt(15.0D0)*0.25D0)
          UNX=DELTAX/DELTA
          UNY=DELTAY/DELTA
          UCTG=dacos(-0.25D0)
        else
          UNZ=0.0D0
          UNX=-dsqrt(0.5D0)
          UNY=dsqrt(0.5D0)
          UCTG=BET
        endif
        do I=1,NPTS
          USOL(1,I,1)=UMX             ! Mx
          USOL(2,I,1)=UMY             ! My
          USOL(3,I,1)=UMZ             ! Mz
          USOL(4,I,1)=UNX             ! Nx   !!
          USOL(5,I,1)=UNY             ! Ny   !!
          USOL(6,I,1)=UNZ             ! Nz
          USOL(7,I,1)=UCTG            ! TETA         !!!!!+/-
        enddo
        do I=1,NPTS
          do J=1,NPDE
            do K=2,3
              USOL(J,I,K)=0D0
            enddo
          enddo
        enddo
        return
      end
C-- USP(X) ----- CSI OF SOLUTION ------------------------------------
      double precision function USP(XI,I)
        include 'vmcw.fh'
        include 'par.fh'
        real*8 USOL, XSOL
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)
        do K=1,NPTS
          USM=dsqrt(USOL(5,K,1)**2+USOL(6,K,1)**2+USOL(4,K,1)**2)
          USOL(4,K,1)=USOL(4,K,1)/USM
          USOL(5,K,1)=USOL(5,K,1)/USM
          USOL(6,K,1)=USOL(6,K,1)/USM
        enddo
        AA=1.0D20
        IK=1
        do II=1,NPTS
          BB=(XSOL(II)-XI)**2
          if(BB.LE.AA)THEN
            AA=BB
            IK=II
          endif
        enddo
        USP=USOL(I,IK,1)
        return
      end
C-- UINIT ------ INITIAL CONDITIONS ---------------------------------
      subroutine UINIT(XI,UI,NPDEI)
        include 'vmcw.fh'
        dimension UI(NPDEI)
        do I=1,NPDEI
          UI(I)=USP(XI,I)
        enddo
        return
      end
C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        include 'vmcw.fh'
        dimension U(NPDE),UX(NPDE),UXX(NPDE),
     *       DFDU(NPDE,NPDE),DFDUX(NPDE,NPDE),DFDUXX(NPDE,NPDE)
        do I=1,NPDE
          do J=1,NPDE
            DFDU(I,J)=0.0D0
            DFDUX(I,J)=0.0D0
            DFDUXX(I,J)=0.0D0
          enddo
        enddo
CCC ATTENTION :   DERIVF IS WRONG
        return
      end
