MODULE nmr

CONTAINS

  !! Computes the NMR lineshape
  !! nu0   -- Larmor frequency, Hz
  !! nub   -- Leggett frequency, Hz
  !! gamma -- intrinsic linewidth, Hz,
  !! fac   -- spectrum margins (in units of gamma)
  !! beta(nbeta) -- array of beta_n in radial coordinates
  !! freq(nspec), spec(nspec) -- arrays for result
  subroutine response(nu0,nub,gamma,fac, &
      nbeta, beta, nspec,freq,spec)

    implicit none
    integer :: nbeta, nspec
    real*8 :: spec(nspec), freq(nspec), beta(nbeta)
    real*8 :: nu0,nub,gamma,fac
    integer :: i,j
    real*8 :: numin,numax,help,nur
    real*8 :: bp,bm,rp,rm,dx,fu
    ! Gauusian quadrature coordinates:
    real*8 :: sp = (3D0 + dsqrt(3D0))/6D0
    real*8 :: sm = (3D0 - dsqrt(3D0))/6D0
    real*8 :: pi = 3.1415926535897932D0

    dx = 1D0/(nbeta-1D0)
    numin = nu0 - fac*gamma
    numax = dsqrt(nu0**2 + nub**2) + fac*gamma
    help = 0.5D0 * (nu0**2 + nub**2)

    do i=1,nspec
       freq(i) = numin + (numax-numin)*dble(i-1)/dble(nspec-1)
       spec(i) = 0D0
       do j=1,nbeta-1
          rp = (dble(j-1)+sp)*dx
          rm = (dble(j-1)+sm)*dx
          bp = sp*beta(j+1)+sm*beta(j)
          bm = sm*beta(j+1)+sp*beta(j)

          nur=SQRT(help+SQRT(help**2-(nu0*nub*COS(bp))**2))
          fu=gamma/(pi*(gamma**2+(freq(i)-nur)**2))
          spec(i) = spec(i) + dx*rp*fu

          nur=SQRT(help+SQRT(help**2-(nu0*nub*COS(bm))**2))
          fu=gamma/(pi*(gamma**2+(freq(i)-nur)**2))
          spec(i) = spec(i) + dx*rm*fu
       enddo
    enddo
  end

END MODULE
