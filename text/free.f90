MODULE energies

  USE general

  USE text

  USE glob

  IMPLICIT NONE

  REAL (KIND=dp), SAVE :: dx
  REAL (KIND=dp), SAVE :: nub,nu0
  REAL (KIND=dp), DIMENSION(0:maxnpt), SAVE :: apsi

  REAL (KIND=dp) :: s3 = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s5 = sqrt(5._dp) ! tmp

  REAL (KIND=dp) chia, vd, xir, de, dar, lsg

  CONTAINS

    subroutine en_surf(a,b,E,Ea,Eb)
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' at surface
      !! parameters used: dar,xir,lsg
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: a,b,E,Ea,Eb
      REAL (KIND=dp) :: nr,nf,nz
      REAL (KIND=dp) :: sin_a, sin_b, cos_a, cos_b, sin2b,cos2b
      E=0
      Ea=0
      Eb=0
      cos_a = cos(a)
      sin_a = sin(a)
      cos_b = cos(b)
      sin_b = sin(b)
      sin2b = sin(2*b)
      cos2b = cos(2*b)

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b

      E = E - 5*dar*(s5*nz*nr-s3*nf)**2/16
      Eb = Eb + 5*dar*(s5*nz*nr-s3*nf)*(s5*cos2b*cos_a+s3*cos_b*sin_a)/8
      Ea = Ea - 5*dar*(s5*nz*nr-s3*nf)*(s5*nz*nf + s3*nr)/8

      ! from bending free energy
      E = E + 4*(2+de)*xir**2*sin_b**2/13
      Eb = Eb + 4*(2+de)*xir**2*sin2b/13

      E = E - 2*lsg*xir**2*sin_b**2/13
      Eb = Eb - 2*lsg*xir**2*sin2b/13
    end subroutine

    subroutine en_bulk(r,a,b,da,db, apsi, vz,vr,vf, lz,lr,lf, w, E,Ea,Eb,Eda,Edb)
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' in the bulk
      !! parameters used: nub/nu0,chia,lo,vd,de,xir
      REAL (KIND=dp) :: r,a,b,da,db,E,Ea,Eb,Eda,Edb
      REAL (KIND=dp) :: apsi, vz,vr,vf, lz,lr,lf, w
      REAL (KIND=dp) :: nz,nr,nf, rzz,rzr,rzf
      REAL (KIND=dp) :: sin_a, sin_b, cos_a, cos_b, cos2b, sin2b
      REAL (KIND=dp) :: con1, con2, help, c,s

      cos_a = cos(a)
      sin_a = sin(a)
      cos_b = cos(b)
      sin_b = sin(b)
      cos2b = cos(2*b)
      sin2b = sin(2*b)

      c=-0.25_dp
      s=SQRT(15.)/4.0_dp

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b
      rzr=(1-c)*nz*nr-s*nf
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2

      E = 0
      Ea = 0
      Eb = 0
      Eda = 0
      Edb = 0

      ! magnetic free energy
      E = E + sin_b**2
      Eb = Eb + sin2b

      ! spin-orbit free energy
      E = E + chia*(nub/nu0 * apsi * sin_b)**2
      Eb = Eb + chia*sin2b*(nub/nu0 * apsi)**2

      ! flow free energy
      E = E - 2*(rzr*vr+rzf*vf+rzz*vz)**2/(5*vd**2)

      help = vr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
           + vf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
           + vz*(-(1-c)*sin2b)
      Eb = Eb - 4*(rzr*vr+rzf*vf+rzz*vz)*help/(5*vd**2)

      help = vr*((1-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
           + vf*((1-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
      Ea = Ea - 4*(rzr*vr+rzf*vf+rzz*vz)*help/(5*vd**2)

      ! vortex free energy
      E = E + lo*w*(rzr*lr+rzf*lf+rzz*lz)**2/5

      help = lr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
           + lf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
           + lz*(-(1-c)*sin2b)
      Eb = Eb + lo*w*4*(rzr*lr+rzf*lf+rzz*lz)*help/10
      help = lr*((1-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
           + lf*((1-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
      Ea = Ea + lo*w*4*(rzr*lr+rzf*lf+rzz*lz)*help/10

      ! bending free energy
      con1 = 4*(4+de)*xir**2/13

      E = E + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2) ! (\nabla n)^2 (?)
      Eda = Eda + con1*2*da*sin_b**2
      Edb = Edb + con1*2*db
      Eb = Eb + con1 * 2*sin_b*cos_b*(da**2 + 1/r**2)

      con2 = -(2+de)*xir**2/26
      help=(s5*sin_a-s3*cos_b*cos_a)*db + &
           (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
           (s5*cos_b*sin_a-s3*cos_a)*sin_b/r  ! s3 \div n + s5 n \rot n (?)

      E = E + con2 * help**2

      Eda = Eda + 2*con2*help * (s5*cos_b*cos_a + s3*sin_a)*sin_b
      Edb = Edb + 2*con2*help * (s5*sin_a - s3*cos_b*cos_a)

      Ea = Ea + 2*con2*help* &
        ( (s5*cos_a + s3*cos_b*sin_a)*db &
        - (s5*cos_b*sin_a - s3*cos_a)*sin_b*da &
        + (s5*cos_b*cos_a + s3*sin_a)*sin_b/r)

      Eb = Eb + 2*con2*help* &
        ( (s3*db - s5*sin_b*da)*sin_b*cos_a &
        + (s5*cos_b*cos_a + s3*sin_a)*cos_b*da &
        + (s5*cos_b*sin_a - s3*cos_a)*cos_b/r &
        - s5*sin_b*sin_a*sin_b/r)
    end subroutine


    subroutine egrad(alpha,beta,e,ga,gb)
      !!! Calculate free energy E and derivatives dE/da(i), dE/db(i)
      !
      ! E = Int e(r, a(r),b(r), da/dr, db/dr,...) r dr
      ! Gaussian quadrature is used for calculating integral. For each
      ! mesh point i = 0:N-1 energies em,ep is calculated in points
      ! i+sm, i+sp (sm/sp = 1 +/- 1/sqrt(3)), and
      ! (em*rm + ep*rp)*dr/2 is added to the sum
      !
      ! Derivatives dE/da(i), dE/db(i) are also needed.
      ! Change DA of a(i) (or b(i)) affects 4 terms in E intergral:
      !   at i-sp, i-sm (for i!=0), i+sp, i+sm (for i!=nmax)
      ! Changes of a in these points are DA*sm, DA*sp, DA*sp, DA*sm
      ! Changes of a' in these points are DA/dr, DA/dr, -DA/dr, -DA/dr
      ! We need to calculate dE/DA = Sum(dE/da * da + dE/da' * da')/DA
      ! We also need r*dr/2 factor as in energy calculation
      !
      ! Strightforward approach is to calculate sum for these 4 points
      !   (Ea*sm*dr - Eda)*r/2 for i+sp, i!=nmax
      !   (Ea*sp*dr - Eda)*r/2 for i+sm, i!=nmax
      !   (Ea*sp*dr + Eda)*r/2 for i-sm, i!=0
      !   (Ea*sm*dr + Eda)*r/2 for i-sp, i!=0
      ! but we can calculate E* only in two points instead of 4
      ! and add some terms to both dE/da(i) and dE/da(i+1)

      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      REAL (KIND=dp) :: rp,rm,bp,bm,ap,am,apsip,apsim,e
      REAL (KIND=dp) :: vzp,vrp,vfp,lzp,lrp,lfp,wp
      REAL (KIND=dp) :: vzm,vrm,vfm,lzm,lrm,lfm,wm
      REAL (KIND=dp) :: da,db
      REAL (KIND=dp) E0,Ea,Eb,Eda,Edb

      REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp
      REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp

      do i=0,nmax
         ga(i)=0.0_dp
         gb(i)=0.0_dp
      enddo
      e=0

      do i=0,nmax-1
         rp=(i+sp)*dx
         rm=(i+sm)*dx
         ap=sp*alpha(i+1)+sm*alpha(i)
         am=sm*alpha(i+1)+sp*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)

         apsip=sp*apsi(i+1)+sm*apsi(i)
         apsim=sm*apsi(i+1)+sp*apsi(i)

         vzp=sp*evz(i+1)+sm*evz(i)
         vzm=sm*evz(i+1)+sp*evz(i)
         vrp=sp*evr(i+1)+sm*evr(i)
         vrm=sm*evr(i+1)+sp*evr(i)
         vfp=sp*evf(i+1)+sm*evf(i)
         vfm=sm*evf(i+1)+sp*evf(i)

         lzp=sp*elz(i+1)+sm*elz(i)
         lzm=sm*elz(i+1)+sp*elz(i)
         lrp=sp*elr(i+1)+sm*elr(i)
         lrm=sm*elr(i+1)+sp*elr(i)
         lfp=sp*elf(i+1)+sm*elf(i)
         lfm=sm*elf(i+1)+sp*elf(i)
         wp=sp*ew(i+1)+sm*ew(i)
         wm=sm*ew(i+1)+sp*ew(i)

         da=(alpha(i+1)-alpha(i))/dx
         db=(beta(i+1)-beta(i))/dx

         call en_bulk(rp, ap,bp,da,db, apsip, vzp,vrp,vfp, &
                  lzp,lrp,lfp, wp, E0,Ea,Eb,Eda,Edb)
         ga(i) = ga(i) + (Ea*sm*dx - Eda)*rp/2.0
         gb(i) = gb(i) + (Eb*sm*dx - Edb)*rp/2.0
         ga(i+1) = ga(i+1) + (Ea*sp*dx + Eda)*rp/2.0
         gb(i+1) = gb(i+1) + (Eb*sp*dx + Edb)*rp/2.0
         e = e + rp*e0*0.5*dx

         call en_bulk(rm, am,bm,da,db, apsim, vzm,vrm,vfm, &
                  lzm,lrm,lfm, wm, E0,Ea,Eb,Eda,Edb)
         ga(i) = ga(i) + (Ea*sp*dx - Eda)*rm/2.0
         gb(i) = gb(i) + (Eb*sp*dx - Edb)*rm/2.0
         ga(i+1) = ga(i+1) + (Ea*sm*dx + Eda)*rm/2.0
         gb(i+1) = gb(i+1) + (Eb*sm*dx + Edb)*rm/2.0
         e = e + rm*E0*0.5*dx

      enddo
      ! surface terms
      call en_surf(alpha(nmax),beta(nmax),E0,Ea,Eb)
      e = e + E0
      ga(nmax) = ga(nmax) + Ea
      gb(nmax) = gb(nmax) + Eb
    end subroutine

    subroutine sfun(n,x,f,g)
      !! wrapper for egrad function for using in the TN
      !! calculate f and g from x values
      !! x as array of both alpha and beta values
      !! g is array of both ga, gb
      IMPLICIT NONE
      INTEGER :: i,n
      REAL (KIND=dp) :: f
      REAL (KIND=dp), DIMENSION(n) :: x,g
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      DO i=1,nmax
         alpha(i)=x(i+1)
         beta(i)=x(i+nmax+1)
      END DO
      alpha(0)=x(1)
      beta(0)=0._dp

      chia = fchia(t,p)
      vd   = fvd(t,p)
      xir  = fxih(t,p,h)/r
      de   = fdelta(t,p)
      dar  = fdar(t,p,r)
      lsg  = 3._dp ! see Fig. 1 in Erkki's paper

      CALL egrad(alpha,beta,f,ga,gb)
      DO i=1,nmax
         g(i+1)=ga(i)
         g(i+nmax+1)=gb(i)
      END DO
      g(1)=ga(0)
    end subroutine


END MODULE energies

