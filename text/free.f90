MODULE energies

  USE general

  USE text

  USE glob

  IMPLICIT NONE

  REAL (KIND=dp), SAVE :: dx
  REAL (KIND=dp), SAVE :: nub,nu0
  REAL (KIND=dp), DIMENSION(0:maxnpt), SAVE :: apsi

  REAL (KIND=dp) :: lsg=3._dp ! see Fig. 1 in Erkki's paper

  REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp ! for Gaussian quadrature
  REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp
  REAL (KIND=dp) :: sq = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s3 = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s5 = sqrt(5._dp) ! tmp

  REAL chia, vd, xir, de, dar

  CONTAINS

    SUBROUTINE ab2n(alpha,beta,nz,nr,nf)
      REAL (KIND=dp) :: alpha,beta,nr,nf,nz
      nr=-SIN(beta)*COS(alpha)
      nf=SIN(beta)*SIN(alpha)
      nz=COS(beta)
    END

    SUBROUTINE sfun(n,x,f,g)
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

      chia=fchia(t,p)
      vd=fvd(t,p)
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      dar=fdar(t,p,r)

      CALL egrad(alpha,beta,f,ga,gb)
      CALL egrad_old(alpha,beta,ga,gb)
      DO i=1,nmax
         g(i+1)=ga(i)
         g(i+nmax+1)=gb(i)
      END DO
      g(1)=ga(0)
    END SUBROUTINE sfun


    FUNCTION esurf(alpha,beta) RESULT(e)
      ! Calculates the surface free energy
      ! E = -5/16 (d/aR) (sqrt5 nz*nr - sqrt3 nf)^2 at i=nmax
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: alpha,beta,nr,nf,nz,e
      call ab2n(alpha, beta, nz, nr, nf)
      e=-5*dar*(s5*nz*nr-s3*nf)**2/16

      ! from bending free energy
      e = e + 4*(2+de)*xir**2*SIN(beta)**2/13
      e = e - 2*lsg*xir**2*SIN(beta)**2/13
    END FUNCTION esurf

    SUBROUTINE egrad_old(alpha,beta,ga,gb)
      ! Calculates the first-order derivatives
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      REAL (KIND=dp) :: nr,nf,nz,help,bn,an,con,bi,bip,bim,rp,rm
      REAL (KIND=dp) :: dap,dam,dbp,dbm,bp,bm,chia
      REAL (KIND=dp) :: db,da,aim,ai,aip,ap,am
      REAL (KIND=dp) :: vd,rzr,rzf,rzz,s,c

      bn=beta(nmax)
      an=alpha(nmax)
      call ab2n(an, bn, nz, nr, nf)
      help=s5*COS(2*bn)*COS(an)+s3*COS(bn)*SIN(an)
      gb(nmax)=gb(nmax)+5*dar*(s5*nz*nr-s3*nf)*help/8
      help=s5*nz*nf+s3*nr
      ga(nmax)=ga(nmax)-5*dar*(s5*nz*nr-s3*nf)*help/8
      gb(nmax)=gb(nmax)+4*(2+de)*xir**2*SIN(2*bn)/13
      gb(nmax)=gb(nmax)-2*lsg*xir**2*SIN(2*bn)/13
!
      gb(0)=0._dp
    END SUBROUTINE egrad_old

    subroutine egr(r,a,b,da,db, apsi, vz,vr,vf, lz,lr,lf, w, E,Ea,Eb,Eda,Edb)
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' at some point

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

      E = E + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2)
      Eda = Eda + con1*2*da*sin_b**2
      Edb = Edb + con1*2*db
      Eb = Eb + con1 * 2*sin_b*cos_b*(da**2 + 1/r**2)

      con2 = -(2+de)*xir**2/26
      help=(s5*sin_a-s3*cos_b*cos_a)*db + &
           (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
           (s5*cos_b*sin_a-s3*cos_a)*sin_b/r

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
    end


    subroutine egrad(alpha,beta,e,ga,gb)
      ! Calculate dE/da(i), dE/db(i)
      ! Change DA of a(i) (or b(i)) affects only 4 terms in E intergral:
      !   at i-sp, i-sm (for i!=0), i+sp, i+sm (for i!=maxn)
      ! Changes of a in these points are DA*sm, DA*sp, DA*sp, DA*sm
      ! Changes of a' in these points are DA/dr, DA/dr, -DA/dr, -DA/dr
      ! We need to calculate dE/DA = Sum(dE/da * da + dE/da' * da')/DA
      ! We also need r*dx/2 factor as in energy calculation
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      REAL (KIND=dp) :: rp,rm,bp,bm,ap,am,apsip,apsim,ep,em,e
      REAL (KIND=dp) :: vzp,vrp,vfp,lzp,lrp,lfp,wp
      REAL (KIND=dp) :: vzm,vrm,vfm,lzm,lrm,lfm,wm
      REAL (KIND=dp) :: da,db
      REAL (KIND=dp) Ea,Eb,Eda,Edb

      e = esurf(alpha(nmax),beta(nmax))

      do i=0,nmax
         ga(i)=0.0_dp
         gb(i)=0.0_dp
         if (i.ne.nmax) then
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

           call egr(rp, ap,bp,da,db, apsip, vzp,vrp,vfp, &
                    lzp,lrp,lfp, wp, ep,Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sm*dx - Eda)*rp/2.0
           gb(i) = gb(i) + (Eb*sm*dx - Edb)*rp/2.0
           call egr(rm, am,bm,da,db, apsim, vzm,vrm,vfm, &
                    lzm,lrm,lfm, wm, em,Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sp*dx - Eda)*rm/2.0
           gb(i) = gb(i) + (Eb*sp*dx - Edb)*rm/2.0

           ! note that in radial coord system we need to sum r*e
           e = e + (rp*ep + rm*em)*0.5*dx

         endif
         if (i.ne.0) then
           rp=(i-sm)*dx
           rm=(i-sp)*dx
           ap=sp*alpha(i)+sm*alpha(i-1)
           am=sm*alpha(i)+sp*alpha(i-1)
           bp=sp*beta(i)+sm*beta(i-1)
           bm=sm*beta(i)+sp*beta(i-1)

           apsip=sp*apsi(i)+sm*apsi(i-1)
           apsim=sm*apsi(i)+sp*apsi(i-1)

           vzp=sp*evz(i)+sm*evz(i-1)
           vzm=sm*evz(i)+sp*evz(i-1)
           vrp=sp*evr(i)+sm*evr(i-1)
           vrm=sm*evr(i)+sp*evr(i-1)
           vfp=sp*evf(i)+sm*evf(i-1)
           vfm=sm*evf(i)+sp*evf(i-1)

           lzp=sp*elz(i)+sm*elz(i-1)
           lzm=sm*elz(i)+sp*elz(i-1)
           lrp=sp*elr(i)+sm*elr(i-1)
           lrm=sm*elr(i)+sp*elr(i-1)
           lfp=sp*elf(i)+sm*elf(i-1)
           lfm=sm*elf(i)+sp*elf(i-1)
           wp=sp*ew(i)+sm*ew(i-1)
           wm=sm*ew(i)+sp*ew(i-1)

           da=(alpha(i)-alpha(i-1))/dx
           db=(beta(i)-beta(i-1))/dx

           call egr(rp, ap,bp,da,db, apsip, vzp,vrp,vfp, &
                    lzp,lrp,lfp, wp, ep,Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sp*dx + Eda)*rp/2.0
           gb(i) = gb(i) + (Eb*sp*dx + Edb)*rp/2.0
           call egr(rm, am,bm,da,db, apsim, vzm,vrm,vfm, &
                    lzm,lrm,lfm, wm, em,Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sm*dx + Eda)*rm/2.0
           gb(i) = gb(i) + (Eb*sm*dx + Edb)*rm/2.0
         endif
      enddo
    end

END MODULE energies

