subroutine calctexture(npttext,textpar,nptspec,specpar,initype, &
     textur,resspec,msglev,apsipar)
  USE general
  USE modu
  USE free
  USE nmr
  USE profiles
  IMPLICIT NONE
  include '../lib/he3.f90h'
  INTEGER :: npttext, nptspec, msglev, iflag
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

  INTEGER, PARAMETER :: maxnpar=2*maxnpt+1
#if USEBTN == 1
  INTEGER, PARAMETER :: nprocs=8
  INTEGER, PARAMETER :: lw=3*maxnpt+3*nprocs + 4*nprocs*nprocs + 7 *(maxnpt*nprocs)
#else
  INTEGER, PARAMETER :: lw=14*maxnpar
#endif
  INTEGER :: n
  REAL (KIND=dp), DIMENSION(0:maxnpt) :: alpha,beta,ga,gb
  REAL (KIND=dp) :: eps,e2,e1,omega,gamma,nu,fac,hei
  REAL (KIND=dp) :: rr,rv,ov,ri,rc,kr
  REAL (KIND=dp) :: f, maxbeta
  REAL (KIND=dp), DIMENSION(maxnpar) :: x,g
  REAL (KIND=dp), DIMENSION(0:nptspec) :: spec, freq
  REAL (KIND=dp), DIMENSION(lw) :: w
  REAL (KIND=dp) c,s,nr,nz,nf, rzr,rzf,rzz


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
  lhv = textpar(8)
  chi=textpar(9)
  nub=textpar(10)

  do i=0,npttext
    apsi(i)=apsipar(i)
  enddo

  if (nptspec > 0) then
    gamma = specpar(1)
    fac = specpar(2)
  endif


  h=2*pi*nu0 * 1000 / he3_gyro ! in Gauss

! set parameters

  ! set lhv
  if (lhv.lt.0D0) then
    lhv = he3_text_lhv(t,p)
  endif

  ! set nub
  if (nub.lt.0D0) then
    nub=he3_nu_b(t,p)/1000
  endif

  ! set lambda/omega
  if (lo.lt.0D0) then
    lo=2.5D0 * he3_text_llh(t,p,omega) / he3_text_a(t,p)
  endif

  chia = chi/(he3_text_a(t,p))
  vd   = dsqrt(0.4D0 * he3_text_a(t,p)/ lhv)
  xir  = he3_text_xih(t,p,h)/r
  de   = he3_text_delta(t,p)
  dar  = he3_text_d(t,p) / (he3_text_a(t,p)*r)
  lsg  = 3._dp ! see Fig. 1 in Erkki's paper


  if (msglev > 0) then
    write (*,*) 'T / Tc',t
    write (*,*) 'pressure / bar',p
    write (*,*) 'Larmor freq. (kHz) =',nu0
    write (*,*) 'Longit. freq. (kHz) =',nub
    write (*,*) 'Field (mT) =',h/10
    write (*,*) 'd/aR =',  dar
    write (*,*) 'xih/R =', xir
    write (*,*) 'delta =', de
    write (*,*) 'vd / Omega R =', vd/(omega*r)
    write (*,*) 'Lambda / Omega =',lo
    write (*,*) 'chi/a =', chia
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
      maxbeta=dacos(1D0/dsqrt(5D0))
      do i=0,nmax
        alpha(i)=pi/3D0
        beta(i)= maxbeta* i/nmax
      enddo
    else
      maxbeta=dacos(-1D0/dsqrt(5D0))
      do i=0,nmax
        alpha(i)=-pi/3D0
        beta(i)= maxbeta* i/nmax
      enddo
    endif
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

#if USEBTN == 1
    call btnez(n,x,f,g, w, lw, sfun, iflag)
#else
    call tn(ierror,n,x,f,g,w,lw,sfun,msglev)
#endif

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

      nr=-dsin(beta(i))*dcos(alpha(i))
      nf=dsin(beta(i))*dsin(alpha(i))
      nz=dcos(beta(i))

      c=-0.25D0 !\cos\theta.
      s=dsqrt(15D0)/4D0 !\sin\theta.
      rzr=(1-c)*nz*nr-s*nf ! H*Rij
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2
    enddo
  endif

  if (nptspec > 0) then ! Calculate NMR lineshape
    call response(beta,nu0,nub,gamma,fac,freq,spec)
    ! return the NMR spectrum
    do i=0,ns
      resspec(i,1) = freq(i)
      resspec(i,2) = spec(i)
    enddo
    if(msglev > 0) then ! Find the highest peak in the spectrum
      call peak(spec,ipos,hei)
      write (*,*) 'Maximum absorption =',hei
      write (*,*) 'Position =',-fac*gamma+ipos*(SQRT(nu0**2+nub**2)-nu0+2*fac*gamma)/ns,' kHz'
    endif
  endif
end subroutine calctexture
