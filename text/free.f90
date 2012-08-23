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
      f=energy_int(alpha,beta)
      CALL egrad(alpha,beta,ga,gb)
      CALL egrad_old(alpha,beta,ga,gb)
      DO i=1,nmax
         g(i+1)=ga(i)
         g(i+nmax+1)=gb(i)
      END DO
      g(1)=ga(0)
    END SUBROUTINE sfun

    !!! Textural free energy at given point
    function energy(r, a,b, apsi, vz,vr,vf, lz,lr,lf, w, da,db) RESULT(e)
      REAL (KIND=dp) r, a,b, apsi, vz,vr,vf, lz,lr,lf, w, da,db, e
      REAL (KIND=dp) nz,nr,nf, rzz,rzr,rzf
      REAL (KIND=dp) sin_a,sin_b,cos_a,cos_b

      REAL (KIND=dp) c,s,help, con1,con2
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp

      sin_a = sin(a)
      sin_b = sin(b)
      cos_a = cos(a)
      cos_b = cos(b)

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b
      rzr=(1-c)*nz*nr-s*nf
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2

      e=0
      ! magnetic free energy
      e = e + sin_b**2

      ! spin-orbit free energy
      e = e + chia*(nub/nu0 * apsi * sin_b)**2

      ! flow free energy
      e = e - 2*(rzr*vr+rzf*vf+rzz*vz)**2/(5*vd**2)

      ! vortex free energy
      e = e + lo*w*(rzr*lr+rzf*lf+rzz*lz)**2/5

      ! bending free energy
      con1 = 4*(4+de)*xir**2/13
      con2 = -(2+de)*xir**2/26

      e = e + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2)

      help=(s5*sin_a-s3*cos_b*cos_a)*db + &
           (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
           (s5*cos_b*sin_a-s3*cos_a)*sin_b/r
      e = e + con2 * help**2
    end

    function energy_int(alpha,beta) RESULT(e)
      ! Integral of textural free energy
      ! Gaussian quadrature is used here and below.
      IMPLICIT NONE
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta
      REAL (KIND=dp) :: rp,rm,bp,bm,ap,am,apsip,apsim,ep,em,e
      REAL (KIND=dp) :: vzp,vrp,vfp,lzp,lrp,lfp,wp
      REAL (KIND=dp) :: vzm,vrm,vfm,lzm,lrm,lfm,wm
      REAL (KIND=dp) :: da,db
      INTEGER :: i

      chia=fchia(t,p)
      vd=fvd(t,p)
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      dar=fdar(t,p,r)

      e=0
      do i=0,nmax-1
         ! we will calculate energy in  i + (3 +/- sqrt(3))/6 points
         rp=(i+sp)*dx
         rm=(i+sm)*dx

         ! interpolate all parameters to these points:
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         ap=sp*alpha(i+1)+sm*alpha(i)
         am=sm*alpha(i+1)+sp*alpha(i)
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

         ! da/dr, db/dr:
         da = (alpha(i+1)-alpha(i))/dx
         db = (beta(i+1)-beta(i))/dx

         ep = energy(rp, ap,bp,apsip,vzp,vrp,vfp,lzp,lrp,lfp,wp, da,db)
         em = energy(rm, am,bm,apsim,vzm,vrm,vfm,lzm,lrm,lfm,wm, da,db)
         ! note that in radial coord system we need to sum r*e
         e = e + (rp*ep + rm*em)*0.5*dx
      enddo
      e=e+esurf(alpha(nmax),beta(nmax))
    end function energy_int


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



    function dga(r,s, a, b, da, db)
      REAL (KIND=dp) dga, r,s, a, b, da, db
      REAL (KIND=dp) sin_a, sin_b, cos_a, cos_b
      cos_a = COS(a)
      sin_a = SIN(a)
      cos_b = COS(b)
      sin_b = SIN(b)
      dga = &
        ( db*(s5*sin_a - s3*cos_a*cos_b) + &
          da*(s5*cos_a*cos_b + s3*sin_a)*sin_b + &
             (s5*cos_b*sin_a - s3*cos_a)*sin_b/r) * &
        ( db*(s5*cos_a + s3*cos_b*sin_a)*r*s - &
             (s5*cos_a*cos_b + s3*sin_a)*sin_b*(r-s) + &
          da*(s3*cos_a - s5*cos_b*sin_a)*sin_b*s*r )
    end

    function dgb(r,s, a, b, da, db)
      REAL (KIND=dp) dgb, r,s, a, b, da, db
      REAL (KIND=dp) sin_a, sin_b, cos_a, cos_b
      cos_a = COS(a)
      sin_a = SIN(a)
      cos_b = COS(b)
      sin_b = SIN(b)
      dgb = &
        ( db*(s5*sin_a - s3*cos_a*cos_b) + &
          da*(s5*cos_a*cos_b + s3*sin_a)*sin_b + &
             (s5*cos_b*sin_a - s3*cos_a)*sin_b/r) * &
        ( da*(s5*cos_a*cos_b + s3*sin_a)*cos_b*r + &
             (s5*cos_b*sin_a - s3*cos_a)*cos_b - &
             (s5*da*sin_b - s3*db)*cos_a*sin_b*r - &
             (s5*sin_a - s3*cos_a*cos_b)*r/s - &
           s5*sin_a*sin_b**2 ) * s
    end

    subroutine egr(r,a,b,da,db, Ea,Eb,Eda,Edb)
      !! Calculate dE/da, dE/db, dE/da', dE/db'

      REAL (KIND=dp) :: r,a,b,da,db,Ea,Eb,Eda,Edb
      REAL (KIND=dp) :: con

      Ea = 0
      Eb = 0
      Eda = 0
      Edb = 0

      ! magnetic free energy (e = sin_b**2)
      Eb = Eb + sin(2*b)

      ! spin-orbit free energy
!      e = e + chia*(nub/nu0 * apsi * sin_b)**2

      ! flow free energy
!      e = e - 2*(rzr*vr+rzf*vf+rzz*vz)**2/(5*vd**2)

      ! vortex free energy
!      e = e + lo*w*(rzr*lr+rzf*lf+rzz*lz)**2/5

      ! bending free energy
!      con1 = 4*(4+de)*xir**2/13
!      con2 = -(2+de)*xir**2/26

!      e = e + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2)

      con = 4*(4+de)*xir**2/13
      Edb = Edb + con*2*db
      Eda = Eda + con*2*da*sin(b)**2
      Eb = Eb + con * sin(2*b)*(da**2 + 1/r**2)

      con=-(2+de)*xir**2/26

!      help=(s5*sin_a-s3*cos_b*cos_a)*db + &
!           (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
!           (s5*cos_b*sin_a-s3*cos_a)*sin_b/r
!      e = e + con2 * help**2

    end


    subroutine egrad(alpha,beta,ga,gb)
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
      REAL (KIND=dp) rp,rm, bp,bm, ap,am, da,db
      REAL (KIND=dp) Ea,Eb,Eda,Edb
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
           da=(alpha(i+1)-alpha(i))/dx
           db=(beta(i+1)-beta(i))/dx
           call egr(rp, ap,bp,da,db, Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sm*dx - Eda)*rp/2.0
           gb(i) = gb(i) + (Eb*sm*dx - Edb)*rp/2.0
           call egr(rm, am,bm,da,db, Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sp*dx - Eda)*rm/2.0
           gb(i) = gb(i) + (Eb*sp*dx - Edb)*rm/2.0
         endif
         if (i.ne.0) then
           rp=(i-sm)*dx
           rm=(i-sp)*dx
           ap=sp*alpha(i)+sm*alpha(i-1)
           am=sm*alpha(i)+sp*alpha(i-1)
           bp=sp*beta(i)+sm*beta(i-1)
           bm=sm*beta(i)+sp*beta(i-1)
           da=(alpha(i)-alpha(i-1))/dx
           db=(beta(i)-beta(i-1))/dx
           call egr(rp, ap,bp,da,db, Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sp*dx + Eda)*rp/2.0
           gb(i) = gb(i) + (Eb*sp*dx + Edb)*rp/2.0
           call egr(rm, am,bm,da,db, Ea,Eb,Eda,Edb)
           ga(i) = ga(i) + (Ea*sm*dx + Eda)*rm/2.0
           gb(i) = gb(i) + (Eb*sm*dx + Edb)*rm/2.0
         endif
      enddo
    end


    SUBROUTINE egrad_old(alpha,beta,ga,gb)
      ! Calculates the first-order derivatives
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      REAL (KIND=dp) :: nr,nf,nz,help,bn,an,con,bi,bip,bim,rp,rm
      REAL (KIND=dp) :: dap,dam,dbp,dbm,bp,bm,chia
      REAL (KIND=dp) :: db,da,aim,ai,aip,ap,am
      REAL (KIND=dp) :: vd,rzr,rzf,rzz,s,c
      REAL (KIND=dp) :: vrp,vfp,vzp,vrm,vfm,vzm
      REAL (KIND=dp) :: lrp,lfp,lzp,lrm,lfm,lzm,wm,wp
      REAL (KIND=dp) :: blp,blm,apsip,apsim
      REAL (KIND=dp) :: cos_a, sin_a, cos_b, sin_b

      dar=fdar(t,p,r)
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)

      bn=beta(nmax)
      an=alpha(nmax)

     ! rp=(i+sp)*dx
     ! rm=(i+sm)*dx
     ! bp=sp*beta(i+1)+sm*beta(i)
     ! bm=sm*beta(i+1)+sp*beta(i)
     ! e=e+0.5*dx*(rp*SIN(bp)**2+rm*SIN(bm)**2)

      call ab2n(an, bn, nz, nr, nf)
      help=s5*COS(2*bn)*COS(an)+s3*COS(bn)*SIN(an)
      gb(nmax)=gb(nmax)+5*dar*(s5*nz*nr-s3*nf)*help/8
      help=s5*nz*nf+s3*nr
      ga(nmax)=ga(nmax)-5*dar*(s5*nz*nr-s3*nf)*help/8
!
      gb(nmax)=gb(nmax)+4*(2+de)*xir**2*SIN(2*bn)/13

      con=-(2+de)*xir**2/26
      DO i=0,nmax
         if (i.ne.0) then
           rp=(i-sm)
           rm=(i-sp)
           ap = sp*alpha(i) + sm*alpha(i-1)
           am = sm*alpha(i) + sp*alpha(i-1)
           bp = sp*beta(i) + sm*beta(i-1)
           bm = sm*beta(i) + sp*beta(i-1)
           dam = alpha(i) - alpha(i-1)
           dbm = beta(i) - beta(i-1)

           gb(i)=gb(i) - con*dgb(rp,-sp, ap, bp, dam, dbm)
           ga(i)=ga(i) - con*dga(rp,-sp, ap, bp, dam, dbm)

           gb(i)=gb(i) - con*dgb(rm,-sm, am, bm, dam, dbm)
           ga(i)=ga(i) - con*dga(rm,-sm, am, bm, dam, dbm)
         endif
         if (i.ne.nmax) then
           rp=(i+sp)
           rm=(i+sm)
           ap = sp*alpha(i+1) + sm*alpha(i)
           am = sm*alpha(i+1) + sp*alpha(i)
           bp = sp*beta(i+1) + sm*beta(i)
           bm = sm*beta(i+1) + sp*beta(i)
           dap = alpha(i+1) - alpha(i)
           dbp = beta(i+1) - beta(i)

           gb(i)=gb(i) + con*dgb(rp,sm, ap, bp, dap, dbp)
           ga(i)=ga(i) + con*dga(rp,sm, ap, bp, dap, dbp)

           gb(i)=gb(i) + con*dgb(rm,sp, am, bm, dap, dbp)
           ga(i)=ga(i) + con*dga(rm,sp, am, bm, dap, dbp)
         endif
      enddo

!
      gb(nmax)=gb(nmax)-2*lsg*xir**2*SIN(2*bn)/13
!
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      vd=fvd(t,p)
      DO i=0,nmax-1

         rp=(i+sp)*dx
         ap=sp*alpha(i+1)+sm*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         vzp=sp*evz(i+1)+sm*evz(i)
         vrp=sp*evr(i+1)+sm*evr(i)
         vfp=sp*evf(i+1)+sm*evf(i)

         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+vfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+vzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)-((3-sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         help=vrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+vfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)-((3-sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)

         rm=(i+sm)*dx
         am=sm*alpha(i+1)+sp*alpha(i)
         bm=sm*beta(i+1)+sp*beta(i)
         vzm=sm*evz(i+1)+sp*evz(i)
         vrm=sm*evr(i+1)+sp*evr(i)
         vfm=sm*evf(i+1)+sp*evf(i)

         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+vfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+vzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)-((3+sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
         help=vrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+vfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)-((3+sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
      END DO
      DO i=1,nmax
         rp=(i-1+sp)*dx
         ap=sp*alpha(i)+sm*alpha(i-1)
         bp=sp*beta(i)+sm*beta(i-1)
         vzp=sp*evz(i)+sm*evz(i-1)
         vrp=sp*evr(i)+sm*evr(i-1)
         vfp=sp*evf(i)+sm*evf(i-1)

         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+vfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+vzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)-((3+sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         help=vrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+vfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)-((3+sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)

         rm=(i-1+sm)*dx
         am=sm*alpha(i)+sp*alpha(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         vzm=sm*evz(i)+sp*evz(i-1)
         vrm=sm*evr(i)+sp*evr(i-1)
         vfm=sm*evf(i)+sp*evf(i-1)

         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+vfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+vzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)-((3-sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
         help=vrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+vfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)-((3-sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
      END DO
!
!
      DO i=0,nmax-1
         rp=(i+sp)*dx
         ap=sp*alpha(i+1)+sm*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         lzp=sp*elz(i+1)+sm*elz(i)
         lrp=sp*elr(i+1)+sm*elr(i)
         lfp=sp*elf(i+1)+sm*elf(i)
         wp=sp*ew(i+1)+sm*ew(i)

         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+lfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+lzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)+lo*wp*((3-sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         help=lrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+lfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)+lo*wp*((3-sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10

         rm=(i+sm)*dx
         am=sm*alpha(i+1)+sp*alpha(i)
         bm=sm*beta(i+1)+sp*beta(i)
         lzm=sm*elz(i+1)+sp*elz(i)
         lrm=sm*elr(i+1)+sp*elr(i)
         lfm=sm*elf(i+1)+sp*elf(i)
         wm=sm*ew(i+1)+sp*ew(i)

         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+lfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+lzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)+lo*wm*((3+sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
         help=lrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+lfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)+lo*wm*((3+sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
      END DO
      DO i=1,nmax
         rp=(i-1+sp)*dx
         ap=sp*alpha(i)+sm*alpha(i-1)
         bp=sp*beta(i)+sm*beta(i-1)
         lzp=sp*elz(i)+sm*elz(i-1)
         lrp=sp*elr(i)+sm*elr(i-1)
         lfp=sp*elf(i)+sm*elf(i-1)
         wp=sp*ew(i)+sm*ew(i-1)

         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+lfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+lzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)+lo*wp*((3+sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         help=lrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+lfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)+lo*wp*((3+sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10

         rm=(i-1+sm)*dx
         am=sm*alpha(i)+sp*alpha(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         lzm=sm*elz(i)+sp*elz(i-1)
         lrm=sm*elr(i)+sp*elr(i-1)
         lfm=sm*elf(i)+sp*elf(i-1)
         wm=sm*ew(i)+sp*ew(i-1)

         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+lfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+lzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)+lo*wm*((3-sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
         help=lrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+lfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)+lo*wm*((3-sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
      END DO

      chia=fchia(t,p)
      DO i=1,nmax-1
         rp=(i+sp)*dx
         rm=(i+sm)*dx
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         !blp=ACOS(-1/4+5/4*COS(bp)**2)
         !blm=ACOS(-1/4+5/4*COS(bm)**2)
         apsip=(3+sq)*apsi(i+1)/6+(3-sq)*apsi(i)/6
         apsim=(3-sq)*apsi(i+1)/6+(3+sq)*apsi(i)/6
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rp*(SIN(bp*2)*0.5*(apsip**2)*sm)
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rm*(SIN(bm*2)*0.5*(apsim**2)*sp)

         rp=(i-1+sp)*dx
         rm=(i-1+sm)*dx
         bp=sp*beta(i)+sm*beta(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         !blp=ACOS(-1/4+5/4*COS(bp)**2)
         !blm=ACOS(-1/4+5/4*COS(bm)**2)
         apsip=(3+sq)*apsi(i)/6+(3-sq)*apsi(i-1)/6
         apsim=(3-sq)*apsi(i)/6+(3+sq)*apsi(i-1)/6
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rp*(SIN(bp*2)*0.5*(apsip**2)*sp)
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rm*(SIN(bm*2)*0.5*(apsim**2)*sm)
      END DO

      gb(0)=0._dp
    END SUBROUTINE egrad_old

    ! rp=(i+sp)*dx
    ! rm=(i+sm)*dx
    ! bp=sp*beta(i+1)+sm*beta(i)
    ! bm=sm*beta(i+1)+sp*beta(i)
    ! gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*sm
    ! gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*sp
    ! rp=(i-1+sp)*dx
    ! rm=(i-1+sm)*dx
    ! bp=sp*beta(i)+sm*beta(i-1)
    ! bm=sm*beta(i)+sp*beta(i-1)
    ! gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*sp
    ! gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*sm

END MODULE energies

