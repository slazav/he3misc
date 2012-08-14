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
      f=energy(alpha,beta)
      CALL egrad(alpha,beta,ga,gb)
      DO i=1,nmax
         g(i+1)=ga(i)
         g(i+nmax+1)=gb(i)
      END DO
      g(1)=ga(0)
    END SUBROUTINE sfun


    FUNCTION energy(alpha,beta) RESULT(e)
      ! Calculates the textural free energy
      IMPLICIT NONE
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta
      REAL (KIND=dp) :: e
      e=emagn(beta)
      e=e+spinorbit(beta)
      e=e+esurf(alpha(nmax),beta(nmax))
      e=e+ebend(alpha,beta)
      e=e+eflow(alpha,beta)
      e=e+evortex(alpha,beta)
    END FUNCTION energy


    FUNCTION emagn(beta) RESULT(e)
      ! Calculates the magnetic free energy
      ! Gaussian quadrature is used here and below.
      ! e = int[ sin(b)**2 ]
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: beta
      REAL (KIND=dp) :: rp,rm,bp,bm,e
      DO i=0,nmax-1
         rp=(i+sp)*dx
         rm=(i+sm)*dx
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         e = e + dx*(rp*SIN(bp)**2+rm*SIN(bm)**2)*0.5
      END DO
    END FUNCTION emagn


    FUNCTION spinorbit(beta) RESULT(e)
      !calculates the spin-orbit free energy
      ! e = int[ chia*(nub/nu0 * apsi * sin(b))**2 ]
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: beta
      REAL (KIND=dp) :: rp,rm,bp,bm,e,chia,blp,blm,apsip,apsim
      chia=fchia(t,p)
      e=0.0_dp
      DO i=0,nmax-1
         rp=(i+sp)*dx
         rm=(i+sm)*dx
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
        ! blp=ACOS(-1/4+5/4*(COS(bp))**2)
        ! blm=ACOS(-1/4+5/4*(COS(bm))**2)
         apsip=sp*apsi(i+1)+sm*apsi(i)
         apsim=sm*apsi(i+1)+sp*apsi(i)
         e = e + dx * chia*((nub/nu0)**2) * & 
           ((apsip**2)*rp*SIN(bp)**2 + (apsim**2)*rm*SIN(bm)**2)*0.5
      END DO
    END FUNCTION spinorbit


    FUNCTION esurf(alpha,beta) RESULT(e)
      ! Calculates the surface free energy
      ! e = 
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: alpha,beta,dar,nr,nf,nz,e
      dar=fdar(t,p,r)
      call ab2n(alpha, beta, nz, nr, nf)
      e=-5*dar*(s5*nz*nr-s3*nf)**2/16
    END FUNCTION esurf


    FUNCTION eflow(alpha,beta) RESULT(e)
      ! Calculates the flow free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta
      REAL (KIND=dp) :: s,c,ap,am,rp,rm,bp,bm,e
      REAL (KIND=dp) :: vd,nr,nf,nz,rzr,rzf,rzz
      REAL (KIND=dp) :: vrp,vfp,vzp,vrm,vfm,vzm
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      e=0.0_dp
      vd=fvd(t,p)
      DO i=0,nmax-1
         vrp=evrp(i)
         vfp=evfp(i)
         vzp=evzp(i)
         rp=(i+sp)*dx
         ap=sp*alpha(i+1)+sm*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e-dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)**2/(5*vd**2)
         vrm=evrm(i)
         vfm=evfm(i)
         vzm=evzm(i)
         rm=(i+sm)*dx
         am=sm*alpha(i+1)+sp*alpha(i)
         bm=sm*beta(i+1)+sp*beta(i)
         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e-dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)**2/(5*vd**2)
      END DO
    END FUNCTION eflow


    FUNCTION evortex(alpha,beta) RESULT(e)
      ! Calculates the vortex free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta
      REAL (KIND=dp) :: s,c,ap,am,rp,rm,bp,bm,e
      REAL (KIND=dp) :: wm,wp,nr,nf,nz,rzr,rzf,rzz
      REAL (KIND=dp) :: lrp,lfp,lzp,lrm,lfm,lzm
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      e=0.0_dp
      DO i=0,nmax-1
         lrp=elrp(i)
         lfp=elfp(i)
         lzp=elzp(i)
         wp=ewp(i)
         rp=(i+sp)*dx
         ap=sp*alpha(i+1)+sm*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         call ab2n(ap, bp, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e+dx*lo*rp*wp*(rzr*lrp+rzf*lfp+rzz*lzp)**2/10
         lrm=elrm(i)
         lfm=elfm(i)
         lzm=elzm(i)
         wm=ewm(i)
         rm=(i+sm)*dx
         am=sm*alpha(i+1)+sp*alpha(i)
         bm=sm*beta(i+1)+sp*beta(i)
         call ab2n(am, bm, nz, nr, nf)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e+dx*lo*rm*wm*(rzr*lrm+rzf*lfm+rzz*lzm)**2/10
      END DO
    END FUNCTION evortex


    FUNCTION ebend(alpha,beta) RESULT(e)
      ! Calculates the bending free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta
      REAL (KIND=dp) :: da,db,e,con,help,xir,de
      REAL (KIND=dp) :: ap,am,rp,rm,bp,bm
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      e=0.0_dp
      con=4*(4+de)*xir**2/13
      DO i=0,nmax-1
         da=alpha(i+1)-alpha(i)
         db=beta(i+1)-beta(i)
         e=e+con*(i+0.5)*db**2
         rp=(i+sp)
         rm=(i+sm)
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         e=e+0.5*con*da**2*(rp*SIN(bp)**2+rm*SIN(bm)**2)
         e=e+0.5*con*(SIN(bp)**2/rp+SIN(bm)**2/rm)
      END DO
      con=-(2+de)*xir**2/26
      DO i=0,nmax-1
         da=alpha(i+1)-alpha(i)
         db=beta(i+1)-beta(i)
         rp=(i+sp)
         rm=(i+sm)
         ap=sp*alpha(i+1)+sm*alpha(i)
         am=sm*alpha(i+1)+sp*alpha(i)
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         help=(s5*SIN(ap)-s3*COS(bp)*COS(ap))*db
         help=help+SIN(bp)*(s5*COS(bp)*COS(ap)+s3*SIN(ap))*da
         help=help+SIN(bp)*(s5*COS(bp)*SIN(ap)-s3*COS(ap))/rp
         e=e+0.5*con*rp*help**2
         help=(s5*SIN(am)-s3*COS(bm)*COS(am))*db
         help=help+SIN(bm)*(s5*COS(bm)*COS(am)+s3*SIN(am))*da
         help=help+SIN(bm)*(s5*COS(bm)*SIN(am)-s3*COS(am))/rm
         e=e+0.5*con*rm*help**2
      END DO
!
      e=e+4*(2+de)*xir**2*SIN(beta(nmax))**2/13
!
      e=e-2*lsg*xir**2*SIN(beta(nmax))**2/13
    END FUNCTION ebend


    function dga(r,s,sin_a, sin_b, cos_a, cos_b, da, db)
      REAL (KIND=dp) dga, r,s,sin_a, sin_b, cos_a, cos_b, da, db
      dga = r*( &
        db*(s5*sin_a - s3*cos_a*cos_b) + &
        da*(s5*cos_a*cos_b + s3*sin_a)*sin_b + &
           (s5*cos_b*sin_a - s3*cos_a)*sin_b/r) * &
        (db*(s5*cos_a + s3*cos_b*sin_a)*s - &
            (s5*cos_a*cos_b + s3*sin_a)*sin_b + &
            (s5*cos_a*cos_b + s3*sin_a)*sin_b*s/r + &
         da*(s3*cos_a - s5*cos_b*sin_a)*sin_b*s )
    end

    function dgb(r,s,sin_a, sin_b, cos_a, cos_b, da, db)
      REAL (KIND=dp) dgb, r,s,sin_a, sin_b, cos_a, cos_b, da, db
      dgb = r*( &
        db*(s5*sin_a - s3*cos_a*cos_b) + &
        da*(s5*cos_a*cos_b + s3*sin_a)*sin_b + &
           (s5*cos_b*sin_a - s3*cos_a)*sin_b/r) * &
        ( s3*cos_a*cos_b - s5*sin_a + &
          da*(s5*cos_a*cos_b + s3*sin_a)*cos_b*s + &
             (s5*cos_b*sin_a - s3*cos_a)*cos_b*s/r + &
          (s3*db - s5*da*sin_b)*cos_a*sin_b*s - s5*sin_a*sin_b**2*s/r )
    end



    SUBROUTINE egrad(alpha,beta,ga,gb)
      ! Calculates the first-order derivatives
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      REAL (KIND=dp) :: nr,nf,nz,help,bn,an,con,bi,bip,bim,rp,rm
      REAL (KIND=dp) :: dap,dam,dbp,dbm,dar,xir,de,bp,bm,chia
      REAL (KIND=dp) :: db,da,aim,ai,aip,ap,am
      REAL (KIND=dp) :: vd,rzr,rzf,rzz,s,c
      REAL (KIND=dp) :: vrp,vfp,vzp,vrm,vfm,vzm
      REAL (KIND=dp) :: lrp,lfp,lzp,lrm,lfm,lzm,wm,wp
      REAL (KIND=dp) :: blp,blm,apsip,apsim
      REAL (KIND=dp) :: cos_a, sin_a, cos_b, sin_b

      dar=fdar(t,p,r)
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      DO i=0,nmax
         ga(i)=0.0_dp
         gb(i)=0.0_dp
      END DO
      bn=beta(nmax)
      an=alpha(nmax)

     ! rp=(i+sp)*dx
     ! rm=(i+sm)*dx
     ! bp=sp*beta(i+1)+sm*beta(i)
     ! bm=sm*beta(i+1)+sp*beta(i)
     ! e=e+0.5*dx*(rp*SIN(bp)**2+rm*SIN(bm)**2)

      DO i=1,nmax-1
         rp=(i+sp)*dx
         rm=(i+sm)*dx
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*sm
         gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*sp
         rp=(i-1+sp)*dx
         rm=(i-1+sm)*dx
         bp=sp*beta(i)+sm*beta(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*sp
         gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*sm
      END DO
      i=nmax
      rp=(i-1+sp)*dx
      rm=(i-1+sm)*dx
      bp=sp*beta(i)+sm*beta(i-1)
      bm=sm*beta(i)+sp*beta(i-1)
      gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*sp
      gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*sm
!
      nr=-SIN(bn)*COS(an)
      nf=SIN(bn)*SIN(an)
      nz=COS(bn)
      help=s5*COS(2*bn)*COS(an)+s3*COS(bn)*SIN(an)
      gb(nmax)=gb(nmax)+5*dar*(s5*nz*nr-s3*nf)*help/8
      help=s5*nz*nf+s3*nr
      ga(nmax)=ga(nmax)-5*dar*(s5*nz*nr-s3*nf)*help/8
!
      gb(nmax)=gb(nmax)+4*(2+de)*xir**2*SIN(2*bn)/13
      con=4*(4+de)*xir**2/13
      DO i=1,nmax-1
         dap=alpha(i+1)-alpha(i)
         dam=alpha(i)-alpha(i-1)
         dbp=beta(i+1)-beta(i)
         dbm=beta(i)-beta(i-1)

         gb(i)=gb(i) - 2*(i+0.5)*con*dbp
         gb(i)=gb(i) + 2*(i-0.5)*con*dbm
!
         rp=(i+sp)
         rm=(i+sm)
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         gb(i)=gb(i) + 0.5*con*dap**2*SIN(2*bp)*rp*sm
         gb(i)=gb(i) + 0.5*con*dap**2*SIN(2*bm)*rm*sp
         ga(i)=ga(i) - con*dap*(rp*SIN(bp)**2+rm*SIN(bm)**2)

         rp=(i-1+sp)
         rm=(i-1+sm)
         bp=sp*beta(i)+sm*beta(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         gb(i)=gb(i) + 0.5*con*dam**2*SIN(2*bp)*rp*sp
         gb(i)=gb(i) + 0.5*con*dam**2*SIN(2*bm)*rm*sm
         ga(i)=ga(i) + con*dam*(rp*SIN(bp)**2+rm*SIN(bm)**2)
!
         rp=(i+sp)
         rm=(i+sm)
         bp=sp*beta(i+1)+sm*beta(i)
         bm=sm*beta(i+1)+sp*beta(i)
         gb(i)=gb(i) + 0.5*con*SIN(2*bp)*sm/rp
         gb(i)=gb(i) + 0.5*con*SIN(2*bm)*sp/rm

         rp=(i-1+sp)
         rm=(i-1+sm)
         bp=sp*beta(i)+sm*beta(i-1)
         bm=sm*beta(i)+sp*beta(i-1)
         gb(i)=gb(i)+0.5*con*SIN(2*bp)*sp/rp
         gb(i)=gb(i)+0.5*con*SIN(2*bm)*sm/rm
      END DO

      i=0
      dap=alpha(i+1)-alpha(i)
      rp=(i+sp)
      rm=(i+sm)
      bp=sp*beta(i+1)+sm*beta(i)
      bm=sm*beta(i+1)+sp*beta(i)
      ga(i)=ga(i)-con*dap*(rp*SIN(bp)**2+rm*SIN(bm)**2)

      i=nmax
      dbm=beta(i)-beta(i-1)
      gb(i)=gb(i)+2*(i-0.5)*con*(dbm)
      rp=(i-1+sp)
      rm=(i-1+sm)
      bp=sp*beta(i)+sm*beta(i-1)
      bm=sm*beta(i)+sp*beta(i-1)
      dam=alpha(i)-alpha(i-1)
      gb(i)=gb(i)+0.5*con*dam**2*rp*SIN(2*bp)*sp
      gb(i)=gb(i)+0.5*con*dam**2*rm*SIN(2*bm)*sm
      ga(i)=ga(i)+con*dam*(rp*SIN(bp)**2+rm*SIN(bm)**2)
!
      gb(i)=gb(i)+0.5*con*SIN(2*bp)*sp/rp
      gb(i)=gb(i)+0.5*con*SIN(2*bm)*sm/rm
!
      con=-(2+de)*xir**2/26
      DO i=0,nmax-1
         rp=(i+sp)
         rm=(i+sm)
         bi=beta(i)
         bip=beta(i+1)
         ai=alpha(i)
         aip=alpha(i+1)

         cos_a = COS(sm*ai + sp*aip)
         sin_a = SIN(sm*ai + sp*aip)
         cos_b = COS(sm*bi + sp*bip)
         sin_b = SIN(sm*bi + sp*bip)

         gb(i)=gb(i) + con*dgb(rp,sm,sin_a, sin_b, cos_a, cos_b, aip-ai, bip-bi)
         ga(i)=ga(i) + con*dga(rp,sm,sin_a, sin_b, cos_a, cos_b, aip-ai, bip-bi)

         cos_a = COS(sp*ai + sm*aip)
         sin_a = SIN(sp*ai + sm*aip)
         cos_b = COS(sp*bi + sm*bip)
         sin_b = SIN(sp*bi + sm*bip)

         gb(i)=gb(i) + con*dgb(rm,sp,sin_a, sin_b, cos_a, cos_b, aip-ai, bip-bi)
         ga(i)=ga(i) + con*dga(rm,sp,sin_a, sin_b, cos_a, cos_b, aip-ai, bip-bi)

      END DO
      DO i=1,nmax
         rp=(i-1+sp)
         rm=(i-1+sm)
         bi=beta(i)
         bim=beta(i-1)
         ai=alpha(i)
         aim=alpha(i-1)

         cos_a = COS(sp*ai + sm*aim)
         sin_a = SIN(sp*ai + sm*aim)
         cos_b = COS(sp*bi + sm*bim)
         sin_b = SIN(sp*bi + sm*bim)

         gb(i)=gb(i) - con*dgb(rp,-sp,sin_a, sin_b, cos_a, cos_b, ai-aim, bi-bim)
         ga(i)=ga(i) - con*dga(rp,-sp,sin_a, sin_b, cos_a, cos_b, ai-aim, bi-bim)

         cos_a = COS(sm*ai + sp*aim)
         sin_a = SIN(sm*ai + sp*aim)
         cos_b = COS(sm*bi + sp*bim)
         sin_b = SIN(sm*bi + sp*bim)

         gb(i)=gb(i) - con*dgb(rm,-sm,sin_a, sin_b, cos_a, cos_b, ai-aim, bi-bim)
         ga(i)=ga(i) - con*dga(rm,-sm,sin_a, sin_b, cos_a, cos_b, ai-aim, bi-bim)

      END DO
!
      gb(nmax)=gb(nmax)-2*lsg*xir**2*SIN(2*bn)/13
!
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      vd=fvd(t,p)
      DO i=0,nmax-1
         vrp=evrp(i)
         vfp=evfp(i)
         vzp=evzp(i)
         rp=(i+sp)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=sp*beta(i+1)+sm*beta(i)
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
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
         vrm=evrm(i)
         vfm=evfm(i)
         vzm=evzm(i)
         rm=(i+sm)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=sm*beta(i+1)+sp*beta(i)
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
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
         vrp=evrp(i-1)
         vfp=evfp(i-1)
         vzp=evzp(i-1)
         rp=(i-1+sp)*dx
         ap=(3+sq)*alpha(i)/6+(3-sq)*alpha(i-1)/6
         bp=sp*beta(i)+sm*beta(i-1)
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
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
         vrm=evrm(i-1)
         vfm=evfm(i-1)
         vzm=evzm(i-1)
         rm=(i-1+sm)*dx
         am=(3-sq)*alpha(i)/6+(3+sq)*alpha(i-1)/6
         bm=sm*beta(i)+sp*beta(i-1)
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
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
         lrp=elrp(i)
         lfp=elfp(i)
         lzp=elzp(i)
         wp=ewp(i)
         rp=(i+sp)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=sp*beta(i+1)+sm*beta(i)
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
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
         lrm=elrm(i)
         lfm=elfm(i)
         lzm=elzm(i)
         wm=ewm(i)
         rm=(i+sm)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=sm*beta(i+1)+sp*beta(i)
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
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
         lrp=elrp(i-1)
         lfp=elfp(i-1)
         lzp=elzp(i-1)
         wp=ewp(i-1)
         rp=(i-1+sp)*dx
         ap=(3+sq)*alpha(i)/6+(3-sq)*alpha(i-1)/6
         bp=sp*beta(i)+sm*beta(i-1)
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
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
         lrm=elrm(i-1)
         lfm=elfm(i-1)
         lzm=elzm(i-1)
         wm=ewm(i-1)
         rm=(i-1+sm)*dx
         am=(3-sq)*alpha(i)/6+(3+sq)*alpha(i-1)/6
         bm=sm*beta(i)+sp*beta(i-1)
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
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
    END SUBROUTINE egrad

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

