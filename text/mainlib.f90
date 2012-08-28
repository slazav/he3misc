subroutine calctexture(npttext,textpar,nptspec,specpar,initype, &
     textur,resspec,msglev,apsipar)
  USE general
  USE free
  USE text
  USE nmr
  USE modu
  USE profiles
  IMPLICIT NONE
  INTEGER :: npttext, nptspec, msglev
  REAL (KIND=dp), DIMENSION(10) :: textpar
  ! 1 - temperature / Tc
  ! 2 - pressure, bar
  ! 3 - larmor frequency, kHz
  ! 4 - cylinder radius, cm
  ! 5 - rotation velocity, rad/s
  ! 6 - omega_v, rad/s
  ! 7 - lambda/omega (s/rad) if >=0
  !     use calculated lambda/omega if == -1
  ! 8 - lambga_HV (kg/(m^3 T^2)) if >= 0
  !     use calculated lambga_HV if == -1
  ! 9 - chi (dimensionless, same scale as in the func that is used to calculate it)
  ! 10 - Leggett frequency, kHz (use 0.5bar vlue if -1)
  REAL (KIND=dp), DIMENSION(2) :: specpar
  ! 1 - half-width of NMR line
  ! 2 - margin for automatic region determination
  INTEGER :: initype
  ! 1 - texture w/o 90 deg peak
  ! 2 - texture with 90 deg peak
  ! 3 - use initial configuration from textur
  ! 4 - use initial configuration from textur w/o minimization
  REAL (KIND=dp), DIMENSION(0:npttext) :: apsipar
  ! apsipar is the   A*Psi=sqrt(2*sin(B_µ/2)) -vector of length npttext+1=number_of_discr_interv.+1
  REAL (KIND=dp), DIMENSION(0:npttext,3) :: textur
  ! columns: r, alpha, beta
  REAL (KIND=dp), DIMENSION(0:nptspec,2) :: resspec
  ! columns: f-f0(kHz), absorption

  INTEGER :: i,ierror,ipos,j,jj,kk,nv
  INTEGER, PARAMETER :: maxnpar=2*maxnpt+1,lw=14*maxnpar
  INTEGER :: n
  REAL (KIND=dp), DIMENSION(0:maxnpt) :: alpha,beta,ga,gb
  REAL (KIND=dp) :: eps,e2,e1,omega,gamma,nu,fac,hei
  REAL (KIND=dp) :: rr,rv,ov,ri,rc,kr
  REAL (KIND=dp) :: f, maxbeta
  REAL (KIND=dp), DIMENSION(maxnpar) :: x,g
  REAL (KIND=dp), DIMENSION(0:nptspec) :: spec
  REAL (KIND=dp), DIMENSION(lw) :: w

  if (npttext > maxnpt) then
    textur(0,1) = -1
    return
  endif

  nmax = npttext
  ns = nptspec

  t = textpar(1)
  p = textpar(2)
  nu0 = textpar(3)
  r = textpar(4)
  omega = textpar(5)
  ov = abs(textpar(6))
  lo = textpar(7)
  flhvfix = textpar(8)/1000 ! convert to program units
  chi=textpar(9)
  nub=textpar(10)

  do i=0,npttext
    apsi(i)=apsipar(i)
  enddo

  if (flhvfix*1000 == -1) then
    flhvfix = flhvtheor(t,p)
  endif

  if (nptspec > 0) then
    gamma = specpar(1)
    fac = specpar(2)
  endif

  ! set nub for 0.5bar
  if (nub.lt.0D0) then
    nub=sqrt(14.46/16.8075*(1-t**2)*(44.2121*t**6-64.5411*t**4+16.9909*t**2+16.862)*1000)
  endif

! Juha's code below

  h=2*pi*nu0/20.4 ! in Gauss

  if (lo == -1) then
    rc=xiglf(t,p)*1.0E-5
    ri=SQRT(6.65E-4/(2*pi*omega))
    lo=6.65E-4*(LOG(ri/rc)-0.75)/(2*pi*fvd(t,p)**2)
  endif

  call set_text_pars(t,p,h)

  if (msglev > 0) then
    write (*,*) 'T / Tc',t
    write (*,*) 'pressure / bar',p
    write (*,*) 'Larmor freq. (kHz) =',nu0
    write (*,*) 'Longit. freq. (kHz) =',nub
    write (*,*) 'Field (mT) =',h/10
    write (*,*) 'd/aR =',fdar(t,p,r)
    write (*,*) 'xih/R =',fxih(t,p,h)/r
    write (*,*) 'delta =',fdelta(t,p)
    write (*,*) 'vd / Omega R =',fvd(t,p)/(omega*r)
    write (*,*) 'Lambda / Omega =',lo
    write (*,*) 'chi/a =',fchia(t,p)
  endif

  dx=1._dp/nmax
  n=2*nmax+1

  ! Initial texture
  if (initype >= 3) then ! from textur parameter
    do i = 0,nmax
      alpha(i) = textur(i,2)*pi/180
      beta(i) = textur(i,3)*pi/180
    enddo
  else ! simple guess
    if (initype == 1) then
      maxbeta = ACOS(1._dp/SQRT(5.))
    else
      maxbeta = ACOS(-1._dp/SQRT(5.))
    endif
    do i=0,nmax
      alpha(i)=pi/3
      beta(i)=maxbeta*i/nmax
    enddo
  endif

  ! Do minimization if needed
  if (initype /= 4) then
    ! Pick the appropriate velocity profile
    if (textpar(6) >= 0) then
      call clusterprofile(r,omega,ov)
    else
      call uniformvortcluster(r,omega,ov)
    endif

    !kr=0.5_dp
    !call twistedstate(r,omega,kr)

    ! Minimization routine
    do i=1,nmax
      x(i+1)=alpha(i)
      x(i+nmax+1)=beta(i)
    enddo
    x(1)=alpha(0)
    call tn(ierror,n,x,f,g,w,lw,sfun,msglev)
    do i=1,nmax
      alpha(i)=x(i+1)
      beta(i)=x(i+nmax+1)
    enddo
    alpha(0)=x(1)
    beta(0)=0._dp

    ! Return the texture
    do i=0,nmax
      textur(i,1) = r*i*dx
      textur(i,2) = alpha(i)*180/pi
      textur(i,3) = beta(i)*180/pi
      ! textur(0,2)=fdar(t,p,r)
      ! textur(1,2)=fa(t,p)
      ! textur(2,2)=fchia(t,p)
      ! textur(3,2)=fdelta(t,p)
      ! textur(4,2)=fxih(t,p,h)
      ! textur(5,2)=nub
      ! textur(i,2) = fchia(t,p)*((nub/nu0)**2)
    enddo
  endif

  if (nptspec > 0) then ! Calculate NMR lineshape
    call response(beta,nu0,nub,gamma,fac,spec)
    ! return the NMR spectrum
    do i=0,ns
      resspec(i,1)=-fac*gamma+i*(SQRT(nu0**2+nub**2)-nu0+2*fac*gamma)/ns
      resspec(i,2) = spec(i)
    enddo
    if(msglev > 0) then ! Find the highest peak in the spectrum
      call peak(spec,ipos,hei)
      write (*,*) 'Maximum absorption =',hei
      write (*,*) 'Position =',-fac*gamma+ipos*(SQRT(nu0**2+nub**2)-nu0+2*fac*gamma)/ns,' kHz'
    endif
  endif
end subroutine calctexture
