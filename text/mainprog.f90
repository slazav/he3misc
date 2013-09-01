PROGRAM hydrostatic
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(8)
  INTEGER, PARAMETER :: npt=200, ns=500
  INTEGER :: i, initype
  REAL (KIND=dp), DIMENSION(10) :: textpar
  REAL (KIND=dp), DIMENSION(2) :: specpar
  REAL (KIND=dp), DIMENSION(0:npt,3) :: textur
  REAL (KIND=dp), DIMENSION(0:ns,2) :: spec
  REAL (KIND=dp), DIMENSION(0:npt) :: apsi
  OPEN (10, FILE='texture.dat')
  OPEN (12, FILE='initials.dat')
  OPEN (13, FILE='spec.dat')
  READ (12,*) textpar(1) ! t
  READ (12,*) textpar(2) ! p
  READ (12,*) textpar(3) ! nu0
  READ (12,*) textpar(4) ! r
  READ (12,*) specpar(1) ! gamma
  READ (12,*) specpar(2) ! fac
  READ (12,*) textpar(5) ! omega
  READ (12,*) textpar(6) ! ov
  READ (12,*) textpar(7) ! lo
  READ (12,*) initype
  READ (12,*) textpar(8) !
  READ (12,*) textpar(9) !chi
  READ (12,*) textpar(10) ! nuB
!
  do i=0,npt
    apsi(i)=0D0
  enddo

  call calctexture(npt,textpar,ns,specpar,initype,textur,spec,1,apsi)
! Save the texture
  WRITE (10,*) '# r', 'alpha', 'beta'
  do i=0,npt
    WRITE (10,*) textur(i,1), textur(i,2), textur(i,3)
  enddo
!
! Save the NMR spectrum
  do i=0,ns
    WRITE (13,*) spec(i,1), spec(i,2)
  enddo
!
! Find the highest peak in the spectrum
END PROGRAM hydrostatic
