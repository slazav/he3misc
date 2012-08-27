MODULE spectra

USE general

IMPLICIT NONE

INTEGER, SAVE :: ns=500

CONTAINS

  SUBROUTINE response(beta,nu0,nub,gamma,fac,spec)
    ! Computes the NMR lineshape
    IMPLICIT NONE
    INTEGER :: i,j
    REAL (KIND=dp), DIMENSION(0:ns) :: spec
    REAL (KIND=dp), DIMENSION(0:nmax) :: beta
    REAL (KIND=dp) :: nu0,nub,gamma,nu,fac,numin,numax,help,nur
    REAL (KIND=dp) :: bp,bm,rp,rm,dx,fu
    REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp
    REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp

    dx=1._dp/nmax
    numin=nu0-fac*gamma
    numax=SQRT(nu0**2+nub**2)+fac*gamma
    help=0.5*(nu0**2+nub**2)
    DO i=0,ns
       nu=numin+i*(numax-numin)/ns
       spec(i)=0._dp
       DO j=0,nmax-1
          rp=(j+sp)*dx
          rm=(j+sm)*dx
          bp=sp*beta(j+1)+sm*beta(j)
          bm=sm*beta(j+1)+sp*beta(j)
          nur=SQRT(help+SQRT(help**2-(nu0*nub*COS(bp))**2))
          fu=gamma/(pi*(gamma**2+(nu-nur)**2))
          spec(i) = spec(i) + dx*rp*fu
          nur=SQRT(help+SQRT(help**2-(nu0*nub*COS(bm))**2))
          fu=gamma/(pi*(gamma**2+(nu-nur)**2))
          spec(i) = spec(i) + dx*rm*fu
       ENDDO
    ENDDO
  END SUBROUTINE response



  SUBROUTINE peak(spec,ipos,hei)
    IMPLICIT NONE
    INTEGER :: i,ipos
    REAL (KIND=dp), DIMENSION(0:ns) :: spec
    REAL (KIND=dp) :: numin,numax,hei
    hei=0._dp
    ipos=0
    DO i=0,ns
       IF (spec(i) > hei) THEN
          hei=spec(i)
          ipos=i
       END IF
    END DO
  END SUBROUTINE peak




END MODULE spectra
