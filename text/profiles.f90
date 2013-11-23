MODULE profiles

USE general
USE free

IMPLICIT NONE

CONTAINS

  SUBROUTINE clusterprofile(r,omega,ov)
    ! Vortex-free velocity profile
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,ov,rr
    DO i=0,nmax
       rr=i*dx*r
       ! flow velocity
       evr(i)=0._dp
       evf(i)=0._dp
       evz(i)=0._dp
       ! vortex direction
       elr(i)=0._dp
       elf(i)=0._dp
       elz(i)=1._dp
       ! 
       ew(i)=2*omega
       IF (rr > r*SQRT(ov/omega)) THEN
          evf(i)=omega*rr-ov*r*r/rr
          ew(i)=0._dp
       ENDIF
    ENDDO
  END SUBROUTINE clusterprofile

  SUBROUTINE uniformvortcluster(r,omega,ov)
    ! uniform vortex cluster
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,ov,rr
    DO i=0,nmax
       rr=i*dx*r
       evr(i)=0._dp
       evf(i)=(omega-ov)*rr
       evz(i)=0._dp
       elr(i)=0._dp
       elf(i)=0._dp
       elz(i)=1._dp
       ew(i)=2*ov
    END DO
  END SUBROUTINE uniformvortcluster


  SUBROUTINE twistedstate(r,omega,kr)
    ! Twisted-state velocity profile
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,kr,rr
    DO i=0,nmax
       rr=i*dx*r
       evr(i)=0._dp
       evf(i)=omega*rr-(kr**2/LOG(1+kr**2))*omega*rr/(1+(kr*rr/r)**2)
       evz(i)=-omega*r*(kr/(LOG(1+kr**2)*(1+(kr*rr/r)**2))-1/kr)
       elr(i)=0._dp
       elf(i)=(kr*rr/r)/SQRT(1+(kr*rr/r)**2)
       elz(i)=1._dp/SQRT(1+(kr*rr/r)**2)
       ew(i)=2*omega*(kr**2/LOG(1+kr**2))*(1+(kr*rr/r)**2)**(-1.5)
    END DO
  END SUBROUTINE twistedstate

END MODULE
