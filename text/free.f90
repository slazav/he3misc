MODULE free

  USE modu

  IMPLICIT NONE

  REAL (KIND=dp) :: s3 = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s5 = sqrt(5._dp) ! tmp

  CONTAINS

    subroutine en_surf(a,b,E,Ea,Eb,da2,db2,Eda2,Edb2)
      
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' at surface
      !! parameters used: de, dar,xir,lsg
      IMPLICIT NONE
     ! INTEGER :: i
      REAL (KIND=dp) :: a,b,E,Ea,Eb,Eda2,Edb2,da2,db2
      REAL (KIND=dp) :: nr,nf,nz
      REAL (KIND=dp) :: sin_a, sin_b, cos_a, cos_b, sin2b,cos2b
    
      E=0.d0
      Ea=0.d0
      Eb=0.d0
      Eda2 = 0.d0
      Edb2 = 0.d0
      cos_a = cos(a)
      sin_a = sin(a)
      cos_b = cos(b)
      sin_b = sin(b)
      sin2b = sin(2*b)
      cos2b = cos(2*b)

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b

      E = E - 5.d0*dar*(s5*nz*nr-s3*nf)**2/16.d0
      Eb = Eb + 5.d0*dar*(s5*nz*nr-s3*nf)*(s5*cos2b*cos_a+s3*cos_b*sin_a)/8.d0
      Ea = Ea - 5.d0*dar*(s5*nz*nr-s3*nf)*(s5*nz*nf + s3*nr)/8.d0
         
      !these terms are now implemented in volume integration
     ! if(1==0) then 
     ! ! from bending free energy
     ! E = E + 4*(2+de)*xir**2*sin_b**2/13
     ! Eb = Eb + 4*(2+de)*xir**2*sin2b/13
     ! endif

      E = E - 2*lsg*xir**2*(sin_b*(sin_b+da2)-sqrt(3.d0/5.d0)*db2)/13
      Eb = Eb - 2*lsg*xir**2*(sin2b+cos_b*da2)/13
      
      Eda2 = Eda2 - 2*lsg*xir**2*sin_b/13
      Edb2 = Edb2 - 2*lsg*xir**2*(-sqrt(3.d0/5.d0))/13
           
    end subroutine

    subroutine en_bulk(r,a,b, da1,da2, db1,db2, apsi, vz,vr,vf, & 
                       lz,lr,lf, w, E,Ea,Eb,Eda1,Eda2,Edb1,Edb2)
      
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' in the bulk
      !! parameters used:
      !!   chia*(nub/nu0)^2, for non-zero apsi
      !!   lo for non-zero rotation
      !!   vd for non-zero flow
      !!   de and xir
      REAL (KIND=dp) :: r,a,b,da1,db1,db2,da2,E,Ea,Eb,Eda1,Edb1,Eda2,Edb2
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


      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b

      c=-0.25_dp !\cos\theta 
      s=SQRT(15.)/4.0_dp !\sin\theta 
      rzr=(1-c)*nz*nr-s*nf ! H*Rij
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2

      E = 0.d0
      Ea = 0.d0
      Eb = 0.d0
      Eda1 = 0.d0
      Edb1 = 0.d0
      Eda2 = 0.d0
      Edb2 = 0.d0
      
      ! magnetic free energy F_DH
      E = E + sin_b**2
      Eb = Eb + sin2b
      
      
      ! spin-orbit free energy
      E = E + chia*(nub/nu0 * apsi * sin_b)**2
      Eb = Eb + chia*sin2b*(nub/nu0 * apsi)**2

      ! flow free energy F_HV
      E = E - 2/5 / (vd**2) * (rzr*vr+rzf*vf+rzz*vz)**2

      help = vr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
           + vf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
           + vz*(-(1-c)*sin2b)
      Eb = Eb - 4/5 / (vd**2) * help * (rzr*vr+rzf*vf+rzz*vz)

      help = vr*((1-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
           + vf*((1-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
      Ea = Ea - 4/5 / (vd**2) * help *(rzr*vr+rzf*vf+rzz*vz)

      ! vortex free energy F_LH
      E = E + lo*w/5 * (rzr*lr+rzf*lf+rzz*lz)**2

      help = lr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
           + lf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
           + lz*(-(1-c)*sin2b)
      Eb = Eb + lo*w*2/5 * help * (rzr*lr+rzf*lf+rzz*lz)
      help = lr*((1-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
           + lf*((1-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
      Ea = Ea + lo*w*2/5 * help * (rzr*lr+rzf*lf+rzz*lz)
      
      ! bending free energy F_G
      !See: .../spinwaves/TextLib_Slava1/Derivations/textcalcn_2D.nb
      con1 = 4*(4+de)*xir**2/13

      !r-direction, (\nabla_i n_j)^2
      !note that original 1D code by Kopu integrated -0.5da**2 analytically
      !to a constant: 
      !sin_b**2 instead of (sin_b**2-0.5d0) 
      if(1==1) then
      E = E + con1*(db1**2 + (sin_b**2)*da1**2 + (sin_b**2)/r**2) 
      Eda1 = Eda1 + con1*2*da1*(sin_b**2)
      else
      E = E + con1*(db1**2 + (sin_b**2-0.5d0)*da1**2 + (sin_b**2)/r**2) 
      Eda1 = Eda1 + con1*2*da1*(sin_b**2-0.5d0)
      endif
      Edb1 = Edb1 + con1*2*db1
      Eb = Eb + con1 * 2*sin_b*cos_b*(da1**2 + 1/r**2)
      !f-direction
      E = E + 1/r**2*con1*(2*sin_b**2*da2+sin_b**2*da2**2+db2**2)
      Eda2 = Eda2 +  1/r**2*con1*(2*sin_b**2+2*sin_b**2*da2)
      Edb2 = Edb2 + 1/r**2*con1*2*db2
      Eb = Eb + 1/r**2*con1*2*sin_b*cos_b*(2*da2+da2**2)

      con2 = -(2+de)*xir**2/26
      !r-direction, s3 \div n + s5 n \rot n
      help=(s5*sin_a-s3*cos_b*cos_a)*db1 + &
           (s5*cos_b*cos_a+s3*sin_a)*sin_b*da1 + &
           (s5*cos_b*sin_a-s3*cos_a)*sin_b/r  
      !f-direction
      help=help+1/r*(-s5*cos_a-s3*cos_b*sin_a)*db2 + &
               1/r*(s5*cos_b*sin_a*sin_b-s3*cos_a*sin_b)*da2 
                    
      
      E = E + con2 * help**2
     
      Eda1 = Eda1 + 2*con2*help * (s5*cos_b*cos_a + s3*sin_a)*sin_b
      Eda2 = Eda2 + 2*con2*help * 1/r*(s5*cos_b*sin_a*sin_b-s3*cos_a*sin_b)

      Edb1 = Edb1 + 2*con2*help * (s5*sin_a - s3*cos_b*cos_a)
      Edb2 = Edb2 + 2*con2*help * 1/r*(-s5*cos_a - s3*cos_b*sin_a)

      
      !r
      Ea = Ea + 2*con2*help* &
        ( (s5*cos_a + s3*cos_b*sin_a)*db1 &
        - (s5*cos_b*sin_a - s3*cos_a)*sin_b*da1 &
        + (s5*cos_b*cos_a + s3*sin_a)*sin_b/r)
      !f
      Ea = Ea + 2*con2*help* &
        (   1/r*(s5*sin_a-s3*cos_b*cos_a)*db2  &
        +   1/r*(s5*cos_b*cos_a*sin_b+s3*sin_a*sin_b)*da2)
      !r
      Eb = Eb + 2*con2*help* &
        ( (s3*db1 - s5*sin_b*da1)*sin_b*cos_a &
        + (s5*cos_b*cos_a + s3*sin_a)*cos_b*da1 &
        + (s5*cos_b*sin_a - s3*cos_a)*cos_b/r &
        - s5*sin_b*sin_a*sin_b/r)
      !f
      Eb = Eb + 2*con2*help* &
        ( 1/r*(s3*sin_b*sin_a)*db2  &
        +   1/r*(s5*cos2b*sin_a-s3*cos_a*cos_b)*da2)
      
      
      !Terms not included in the vector expression of bending energy
      !(Hakonen et al 83 eq. (28), second line), they are originally (1D code by Kopu)
      !implemented as surface terms, which angle independent part of them becomes
      !after integration. Constants copyed from original code.
      !For derivation, see .../spinwaves/TextLib_Slava1/Derivations/textcalcn_2D.nb
      E = E + 4.d0*(2+de)*xir**2/13.d0*sin2b*(-db2*da1+(1+da2)*db1)/r
      Eb = Eb + 4.d0*(2+de)*xir**2/13.d0*2*cos2b*(-db2*da1+(1+da2)*db1)/r
      Edb1 = Edb1 + 4.d0*(2+de)*xir**2/13.d0*sin2b*(1+da2)/r
      Edb2 = Edb2 + 4.d0*(2+de)*xir**2/13.d0*sin2b*(-da1)/r
      Eda1 = Eda1 +  4.d0*(2+de)*xir**2/13.d0*sin2b*(-db2)/r
      Eda2 = Eda2 +  4.d0*(2+de)*xir**2/13.d0*sin2b*(da2*db1)/r
         

    end subroutine


    subroutine egrad(n,alpha,beta,e,ga,gb)
      !!! Calculate integral free energy E and derivatives dE/da(i), dE/db(i)
      !
      ! alpha(0) - central point, alpha(i_r+(i_f-1)*nrmax) - (i_r, i_f) point

      !!!!!!COMMENTS BELOW ARE OUTDATED (1d version) !
      ! E = Int e(r, a(r),b(r), da/dr, db/dr,...) r dr
      ! Gaussian quadrature is used for calculating integral. For each
      ! mesh point i = 0:N-1 energies em,ep is calculated in points
      ! i+sm, i+sp (sm/sp = (3 +/- 3/sqrt(3))/6), and
      ! (em*rm + ep*rp)*dr/2 is added to the sum
      !
      ! Derivatives dE/da(i), dE/db(i) are also needed.
      ! Change DA in a(i) (or b(i)) affects 4 terms in E intergral:
      !   at i-sp, i-sm (for i!=0), i+sp, i+sm (for i!=n)
      ! Changes of a in these points are DA*sm, DA*sp, DA*sp, DA*sm
      ! Changes of a' in these points are DA/dr, DA/dr, -DA/dr, -DA/dr
      ! We need to calculate dE/DA = Sum(dE/da * da + dE/da' * da')/DA
      ! We also need r*dr/2 factor as in energy calculation
      !
      ! Strightforward approach is to calculate sum for these 4 points
      !   (Ea*sm*dr - Eda)*r/2 for i+sp, i!=n
      !   (Ea*sp*dr - Eda)*r/2 for i+sm, i!=n
      !   (Ea*sp*dr + Eda)*r/2 for i-sm, i!=0
      !   (Ea*sm*dr + Eda)*r/2 for i-sp, i!=0
      ! but we can calculate E* only in two points instead of 4
      ! and add some terms to both dE/da(i) and dE/da(i+1)
      !!!!!!DOWN TO HERE


      IMPLICIT NONE
      INTEGER :: i,ii,n,k,kk,kkk,in1,in2
      INTEGER, DIMENSION(1:4,1:2) :: ind
      REAL (KIND=dp), DIMENSION(0:nmaxt) :: alpha,beta,ga,gb
      REAL (KIND=dp), DIMENSION(0:nrmax,1:nfmax+1) :: alpha_mat,beta_mat,apsi_mat,ga_mat,gb_mat
      REAL (KIND=dp) :: bpm,apm,apsipm,e,ap,bp,am,bm
      REAL (KIND=dp), DIMENSION(1:2) :: vzpm,vrpm,vfpm,lzpm,lrpm,lfpm,wpm
     
      REAL (KIND=dp) :: da1,db1,da2,db2
      REAL (KIND=dp) :: E0,Ea,Eb,Eda1,Edb1,Eda2,Edb2
      REAL (KIND=dp),DIMENSION(1:2) :: rpm
      REAL (KIND=dp),DIMENSION(1:4,1:2,1:2) :: spm,gr1,gr2

      REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp
      REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp

      REAL (KIND=dp) :: df, nfhelp, testhelp, help1
      nfhelp=1.d0/REAL(nfmax)
      df=2.d0*pi*nfhelp
      
      !periodic boundary conditions->n_intervals=n_points
      !construct help arrays

      !help indices for r and f directions in interpolation:
      !to take steps over corners (real grid points)
      ind(:,1)= (/ 0,1,0,1 /) !r-direction
      ind(:,2)= (/ 0,0,1,1 /) !f-direction
      

      !weights for bilinear interpolation to gaussian quadrature points (4)
      spm(:,1,1) = (/ sp*sp,sp*sm,sp*sm,sm*sm /) !to first interp. point from 4 corner points
      spm(:,2,1) = (/ sp*sm,sp*sp,sm*sm,sp*sm /) !second interp. point
      spm(:,1,2) = (/ sp*sm,sm*sm,sp*sp,sp*sm /) !third  interp. point
      spm(:,2,2) = (/ sm*sm,sm*sp,sm*sp,sp*sp /) !fourth interp. point
      
      !weights for gradients of the interpolation (1 is r direction,2 is f direction)
      gr1(:,1,1) = (/ -sp,sp,-sm,sm /) !first interp. point from 4 corner points
      gr1(:,2,1) = (/ -sp,sp,-sm,sm /) !second interp. point
      gr1(:,1,2) = (/ -sm,sm,-sp,sp /) !third interp. point
      gr1(:,2,2) = (/ -sm,sm,-sp,sp /) !fourth interp. point
      gr1=gr1/dx
     
      gr2(:,1,1) = (/ -sp,-sm,sp,sm /) !first interp. point
      gr2(:,2,1) = (/ -sm,-sp,sm,sp /) !second interp. point
      gr2(:,1,2) = (/ -sp,-sm,sp,sm /) !third interp. point
      gr2(:,2,2) = (/ -sm,-sp,sm,sp /) !fourth interp. point
      gr2=gr2/df


      
      alpha_mat=0.d0
      beta_mat=0.d0
      apsi_mat=0.d0
      ga_mat=0.d0
      gb_mat=0.d0      

      do ii=1,nfmax
         do i=1,nrmax
            alpha_mat(i,ii) = alpha(i+(ii-1)*nrmax)
            beta_mat(i,ii) = beta(i+(ii-1)*nrmax)
            apsi_mat(i,ii) = apsi(i+(ii-1)*nrmax)
         enddo
         alpha_mat(0,ii) = alpha(0)
         beta_mat(0,ii) = beta(0)
         apsi_mat(0,ii) = apsi(0)  
      enddo
      
      !for periodic boundaries in f direction, copy first column to the nfmax+1:th column
      do i=0,nrmax         
         alpha_mat(i,nfmax+1) = alpha_mat(i,1)
         beta_mat(i,nfmax+1) = beta_mat(i,1)
         apsi_mat(i,nfmax+1) = apsi_mat(i,1)         
      enddo
      

     
     
     e=0.d0
  
     Eda1=0.d0
     Eda2=0.d0
     Edb1=0.d0
     Edb2=0.d0

   
      
      do i=0,nrmax-1
         !f-independent iterpolation
         rpm(2)=(i+sp)*dx
         rpm(1)=(i+sm)*dx
                 
         vzpm(2)=sp*evz(i+1)+sm*evz(i)
         vzpm(1)=sm*evz(i+1)+sp*evz(i)
         vrpm(2)=sp*evr(i+1)+sm*evr(i)
         vrpm(1)=sm*evr(i+1)+sp*evr(i)
         vfpm(2)=sp*evf(i+1)+sm*evf(i)
         vfpm(1)=sm*evf(i+1)+sp*evf(i)

         lzpm(2)=sp*elz(i+1)+sm*elz(i)
         lzpm(1)=sm*elz(i+1)+sp*elz(i)
         lrpm(2)=sp*elr(i+1)+sm*elr(i)
         lrpm(1)=sm*elr(i+1)+sp*elr(i)
         lfpm(2)=sp*elf(i+1)+sm*elf(i)
         lfpm(1)=sm*elf(i+1)+sp*elf(i)
         wpm(2)=sp*ew(i+1)+sm*ew(i)
         wpm(1)=sm*ew(i+1)+sp*ew(i)

      !f-dependent interpolation
      do ii=1,nfmax
         
         !periodic boundaries: first point at 0.5*df 
         !(phi is not needed)
         !phip=(ii+sp)*df-0.5d0*df
         !phim=(ii+sm)*df-0.5d0*df


      !k and kk run over interpolation grid: k in r-direction, kk in phi-direction
      do k=1,2
      do kk=1,2
        
	apm=0.d0
	bpm=0.d0
	apsipm=0.d0
	da1=0.d0
	da2=0.d0
	db1=0.d0
	db2=0.d0
      do kkk=1,4 !kkk runs over real grid points (corners)

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)

         !example:
         !apm =spm(1,k,kk)*alpha_mat(i,ii)+spm(2,k,kk)*alpha_mat(i+1,ii) & 
         !     +spm(3,k,kk)*alpha_mat(i,ii+1)+spm(4,k,kk)*alpha_mat(i+1,ii+1)
         
         apm = apm + spm(kkk,k,kk)*alpha_mat(in1,in2)                  
         bpm = bpm + spm(kkk,k,kk)*beta_mat(in1,in2)                       
         apsipm = apsipm + spm(kkk,k,kk)*apsi_mat(in1,in2)
        
         da1 = da1 + (gr1(kkk,k,kk)*alpha_mat(in1,in2))
         db1 = db1 + (gr1(kkk,k,kk)*beta_mat(in1,in2))
         
         da2 = da2 + (gr2(kkk,k,kk)*alpha_mat(in1,in2))
         db2 = db2 + (gr2(kkk,k,kk)*beta_mat(in1,in2))

         
         
      enddo !Loop over 4 grid points (corners)
	!cut loop to include all 4 corner points at each interpolation point: en_bulk is NONLINEAR
     
     

     
      
      Call en_bulk(rpm(k),apm,bpm,da1,da2,db1,db2,apsipm, vzpm(k),vrpm(k),vfpm(k), &
                  lzpm(k),lrpm(k),lfpm(k), wpm(k), E0,Ea,Eb,Eda1,Eda2,Edb1,Edb2)
         !en_bulk call sequence:(r,a,b, da1,da2, db1,db2, apsi, vz,vr,vf, & 
         !        lz,lr,lf, w, E,Ea,Eb,Eda1,Eda2,Edb1,Edb2)
     
     
	!write energy: not iterated over grid points
      e = e + rpm(k)*e0

      

      !help1=0.d0
      do kkk=1,4

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)       
                 
        
   
         ga_mat(in1,in2)=ga_mat(in1,in2) + (Ea*spm(kkk,k,kk) + Eda1*gr1(kkk,k,kk) &
              + Eda2*gr2(kkk,k,kk))*rpm(k)

         gb_mat(in1,in2)=gb_mat(in1,in2) + (Eb*spm(kkk,k,kk) + Edb1*gr1(kkk,k,kk) &
              + Edb2*gr2(kkk,k,kk))*rpm(k)
         
        testhelp = (Ea*spm(kkk,k,kk) + Eda1*gr1(kkk,k,kk)+ Eda2*gr2(kkk,k,kk))*rpm(k)*dx*df*0.25d0  
       
      enddo !loop over 4 grid points (corners)


      enddo !quadrature, f-direction
      enddo !quadrature, r-direction

      enddo !f-loop
      enddo !r-loop
      !set correct scale
      e=e*0.25d0*dx*df
      ga_mat=ga_mat*0.25d0*dx*df
      gb_mat=gb_mat*0.25d0*dx*df    

      
      
      ! surface terms: NOTE THAT RADIUS IS INCLUDED IN en_surf, no need to multiply by r
      !furthermore, Rmax=1 in this code
         
     
      in1=nrmax !boundary
      
      do ii=1,nfmax
        ap=sp*alpha_mat(in1,ii+1)+sm*alpha_mat(in1,ii)	
        am=sm*alpha_mat(in1,ii+1)+sp*alpha_mat(in1,ii)
        bp=sp*beta_mat(in1,ii+1)+sm*beta_mat(in1,ii)
        bm=sm*beta_mat(in1,ii+1)+sp*beta_mat(in1,ii)
	da2=(alpha_mat(in1,ii+1)-alpha_mat(in1,ii))/df
	db2=(beta_mat(in1,ii+1)-beta_mat(in1,ii))/df
      

      call en_surf(ap,bp,E0,Ea,Eb,da2,db2,Eda2,Edb2)
     
      
      !en_surf call sequence:(a,b,E,Ea,Eb,da2,db2,Eda2,Edb2)
      e = e + E0*df*0.5d0 !no need to multiply by radius,
      !radius included in en_surf (and total radius is unity)
      ga_mat(in1,ii) = ga_mat(in1,ii) + (Ea*sm*df - Eda2)*0.5d0
      gb_mat(in1,ii) = gb_mat(in1,ii) + (Eb*sm*df - Edb2)*0.5d0
      ga_mat(in1,ii+1) = ga_mat(in1,ii+1) + (Ea*sp*df + Eda2)*0.5d0
      gb_mat(in1,ii+1) = gb_mat(in1,ii+1) + (Eb*sp*df + Edb2)*0.5d0

      call en_surf(am,bm,E0,Ea,Eb,da2,db2,Eda2,Edb2)
     
      
      !en_surf call sequence:(a,b,E,Ea,Eb,da2,db2,Eda2,Edb2)
      e = e + E0*df*0.5d0 !no need to multiply by radius, 
      !radius included in en_surf (and total radius is unity)
      ga_mat(in1,ii) = ga_mat(in1,ii) + (Ea*sp*df - Eda2)*0.5d0
      gb_mat(in1,ii) = gb_mat(in1,ii) + (Eb*sp*df - Edb2)*0.5d0
      ga_mat(in1,ii+1) = ga_mat(in1,ii+1) + (Ea*sm*df + Eda2)*0.5d0
      gb_mat(in1,ii+1) = gb_mat(in1,ii+1) + (Eb*sm*df + Edb2)*0.5d0

      enddo
      
     
     
               
      ga=0.d0
      gb=0.d0
      
      !copy gradients from help slots
      do i=0,nrmax 
        ga_mat(i,1)=ga_mat(i,1)+ga_mat(i,nfmax+1)
        gb_mat(i,1)=gb_mat(i,1)+gb_mat(i,nfmax+1)                 
      enddo
      !reshape gradient matrices into vectors
      do ii=1,nfmax
      do i=1,nrmax
      ga(i+(ii-1)*nrmax) =  ga_mat(i,ii)
      gb(i+(ii-1)*nrmax) =  gb_mat(i,ii)     
      enddo
      ga(0) = ga(0) + ga_mat(0,ii)
      gb(0) = gb(0) + gb_mat(0,ii)
      enddo    
     
    
    end subroutine

    subroutine sfun(n,x,f,g)
      !! Wrapper for egrad function for using in the TN.
      !! Calculate f and g from x values.
      !! x as array of both alpha and beta values
      !! g is array of both ga, gb
      IMPLICIT NONE
      INTEGER :: i,n !,nmax
      REAL (KIND=dp) :: f, e1, ea1, eb1, eda1, edb1, eda2, edb2
      REAL (KIND=dp), DIMENSION(n) :: x,g
      REAL (KIND=dp), DIMENSION(0:nmaxt) :: alpha,beta,ga,gb
      !nmax=(n-1)/2
      do i=0,nmaxt
         alpha(i)=x(i+1)
         beta(i)=x(i+nmaxt+2)
      enddo
      
      call egrad(nmaxt,alpha,beta,f,ga,gb)      
      !egrad call sequence: (n,alpha,beta,e,ga,gb)
    
      do i=0,nmaxt
         g(i+1)=ga(i)
         g(i+nmaxt+2)=gb(i) 
         !WRITE(25,*) i, alpha(i), beta(i), ga(i), gb(i)
         
      enddo
      
         !WRITE(26,*) 0,0,0,0,0,0,0,0,0,0,0

    end subroutine


END MODULE

