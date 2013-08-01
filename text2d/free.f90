MODULE free

  USE modu
  USE chkder

  IMPLICIT NONE

  REAL (KIND=dp) :: s3 = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s5 = sqrt(5._dp) ! tmp
 

  CONTAINS
    
    !subroutine en_center(b,E,Eb)
    !  REAL (KIND=dp) :: b,E,Eb, sin_b, sin2b
    !  
    !  sin2b = sin(2*b)
    !  sin_b = sin(b)
    !  
    !
    !  E = -2._dp*pi*4._dp*(2._dp+de)*xir**2*sin_b**2/13._dp
    !  Eb =  - 2._dp*pi*4._dp*(2._dp+de)*xir**2*sin2b/13._dp
    !  
    !end subroutine en_center


    subroutine en_surf(a,b,E,Ea,Eb,da2,db2,Eda2,Edb2)
      
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' at surface
      !! parameters used: de, dar,xir,lsg
      IMPLICIT NONE
     ! INTEGER :: i
      REAL (KIND=dp) :: a,b,E,Ea,Eb,Eda2,Edb2,da2,db2
      REAL (KIND=dp) :: nr,nf,nz
      REAL (KIND=dp) :: sin_a, sin_b, cos_a, cos_b, sin2b,cos2b
    
      E=0._dp
      Ea=0._dp
      Eb=0._dp
      Eda2 = 0._dp
      Edb2 = 0._dp
      cos_a = cos(a)
      sin_a = sin(a)
      cos_b = cos(b)
      sin_b = sin(b)
      sin2b = sin(2*b)
      cos2b = cos(2*b)

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b

      E = E - 5._dp*dar*(s5*nz*nr-s3*nf)**2/16._dp
      Eb = Eb + 5._dp*dar*(s5*nz*nr-s3*nf)*(s5*cos2b*cos_a+s3*cos_b*sin_a)/8._dp
      Ea = Ea - 5._dp*dar*(s5*nz*nr-s3*nf)*(s5*nz*nf + s3*nr)/8._dp
         
      
     ! if(1==0) then  !surface contribution, from volume gradient energy
     ! ! from bending free energy
     ! E = E + 4._dp*(2._dp+de)*xir**2*sin_b**2/13._dp
     ! Eb = Eb + 4._dp*(2._dp+de)*xir**2*sin2b/13._dp
     ! else if(1==0) then !testing
     !    E = E + 0.5_dp*xir**2*sin_b**2
     ! endif

      E = E - 2._dp*lsg*xir**2*sin_b*(sin_b*(1-da2)+sqrt(3._dp/5._dp)*db2)/13._dp
      Eb = Eb - 2._dp*lsg*xir**2*cos_b*(2*sin_b*(1-da2)+sqrt(3._dp/5._dp)*db2)/13._dp
      
      Eda2 = Eda2 + 2._dp*lsg*xir**2*sin_b**2/13._dp
      Edb2 = Edb2 - 2._dp*lsg*xir**2*sin_b*sqrt(3._dp/5._dp)/13._dp
           
    end subroutine

    subroutine en_bulk(r,f,nr,nf,nz,dnr1,dnr2,dnf1,dnf2,dnz1,dnz2, apsi, & 
         vz,vr,vf, lz,lr,lf, w, & 
         E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2)      
      
      !! parameters used:
      !!   chia*(nub/nu0)^2, for non-zero apsi
      !!   lo for non-zero rotation
      !!   vd for non-zero flow
      !!   de and xir
      REAL (KIND=dp) :: r,f,E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2
      REAL (KIND=dp) :: apsi, vz,vr,vf, lz,lr,lf, w
      REAL (KIND=dp) :: nz,nr,nf, rzz,rzr,rzf, dnr1,dnr2,dnf1,dnf2,dnz1,dnz2
     
      REAL (KIND=dp) :: con1, con2, help, c,s, cons, cons2, helpr, helpf, helpz, c1, c2
      REAL (KIND=dp) :: hr, hf, hz, hbeta, halpha
      REAL (KIND=dp) :: helpt1, helpt2

      c1 = 5._dp/4._dp ! tmp
      c2 = -sqrt(15._dp)/4._dp ! tmp

     
      c=-0.25_dp !\cos\theta 
      s=SQRT(15.)/4.0_dp !\sin\theta 
      rzr=(1-c)*nz*nr-s*nf ! H*Rij
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2

      E = 0._dp

      Er = 0._dp
      Ef = 0._dp
      Ez = 0._dp
      Edr1 = 0._dp
      Edf1 = 0._dp
      Edr2 = 0._dp
      Edf2 = 0._dp
      Edz1 = 0._dp
      Edz2 = 0._dp


      ! magnetic free energy F_DH

      if(1==0) then !apply tilt
            
      hbeta=pi/180._dp*5._dp
      halpha=0._dp

      hr=-sin(hbeta)*cos(halpha+f)
      hf=sin(hbeta)*sin(halpha+f)
      hz=cos(hbeta)
      
     
      help=hr*nr + hf*nf + hz*nz
     

      else !no tilt
      E = E - (nz**2-1._dp)
      Ez = Ez - 2._dp*nz
      endif
     


      ! spin-orbit free energy
      E = E - chia*(nub/nu0*apsi)**2*(nz**2-1._dp)
      Ez = Ez - 2._dp*chia*nz*(nub/nu0*apsi)**2
      

      
      ! flow free energy F_HV  
      help=(rzr*vr+rzf*vf+rzz*vz)
      cons=- 2._dp/5._dp / (vd**2)
      E = E +cons*help**2
     !rzr=(1-c)*nz*nr-s*nf ! H*Rij
     !rzf=(1-c)*nz*nf+s*nr
     !rzz=c+(1-c)*nz**2
      Er = Er + 2._dp*cons*help*((1-c)*nz*vr+s*vf)
      Ef = Ef + 2._dp*cons*help*((1-c)*nz*vf-s*vr)
      Ez = Ez + 2._dp*cons*help*(1-c)*(2._dp*nz*vz+vr*nr+vf*nf)

      ! vortex free energy F_LH
      help=(rzr*lr+rzf*lf+rzz*lz)
      cons2=lo*w/5._dp
      E = E + cons2* help**2
     
      Er = Er + 2._dp*cons2*help*((1-c)*nz*lr+s*lf)
      Ef = Ef + 2._dp*cons2*help*((1-c)*nz*lf-s*lr)
      Ez = Ez + 2._dp*cons2*help*(1-c)*(2._dp*nz*lz+lr*nr+lf*nf)
      
     
      ! bending free energy F_G
      !See: .../spinwaves/.../Derivations/textcalcn_2D_CORRECT.nb and textcalcn_2D_no_angles.nb
      !Here we use the tensor expression directly (transformed into vector component form)


      !see textcalcn_2D_CORRECT.nb
      !con1 = 2._dp*(4._dp+de)*xir**2/13._dp      
      con1 = 8._dp*xir**2/13._dp    
      ! con1 = 8._dp*(9._dp+2._dp*de)*xir**2/65._dp
      !! con1= 0.5_dp*xir**2
      


      !(\nabla_i n_j)^2, see textcalcn_2D_deln_sq.nb
      E = E + con1*(dnr1**2+dnf1**2+dnz1**2 + & 
           ((nr+dnf2)**2+(nf-dnr2)**2+dnz2**2)/r**2)
      
      Er = Er + 2._dp*con1*(nr+dnf2)/r**2
      Ef = Ef + 2._dp*con1*(nf-dnr2)/r**2

      Edr1 = Edr1 + 2._dp*con1*dnr1
      Edf1 = Edf1 + 2._dp*con1*dnf1
      Edz1 = Edz1 + 2._dp*con1*dnz1
      
      Edr2 = Edr2 - 2._dp*con1*(nf-dnr2)/r**2
      Edf2 = Edf2 + 2._dp*con1*(nr+dnf2)/r**2
      Edz2 = Edz2 + 2._dp*con1*dnz2/r**2
      
      
     
      con2 = 8._dp*(2._dp+de)*xir**2/65._dp  !see textcalcn_2D_CORRECT.nb
      !!con2 = xir**2/5._dp

      ! {5/4*(n /div n +(n /nabla) n)-sqrt(15)/4 /rot(n)}^2}^2 , see "material derivative" for (n /nabla)
      helpr = c1*(nr*(dnr1+dnf2/r+nr/r) + nr*dnr1+nf*dnr2/r-nf**2/r) + c2*dnz2/r
      helpf = c1*(nf*(dnr1+dnf2/r+nr/r) + nr*dnf1+nf*dnf2/r+nr*nf/r) + c2*(-dnz1)
      helpz = c1*(nz*(dnr1+dnf2/r+nr/r) + nr*dnz1+nf*dnz2/r) + c2*(nf+dnf1*r-dnr2)/r

      E = E +con2*(helpr**2+helpf**2+helpz**2)
         

      Er = Er + 2._dp*con2*c1*(helpr*(2._dp*dnr1+dnf2/r+2._dp*nr/r)+helpf*(2._dp*nf/r+dnf1)+helpz*(nz/r+dnz1))
      Ef = Ef + 2._dp*con2*(helpr*c1*(dnr2/r-2._dp*nf/r)+helpf*c1*(dnr1+2._dp*dnf2/r+2._dp*nr/r) +helpz*(c1*dnz2/r+c2/r))
      Ez = Ez + 2._dp*con2*helpz*c1*(dnr1+dnf2/r+nr/r)

      Edr1 = Edr1 + 2._dp*con2*c1*(helpr*2._dp*nr+helpf*nf+helpz*nz)
      Edf1 = Edf1 + 2._dp*con2*(helpz*c2+c1*helpf*nr)
      Edz1 = Edz1 + 2._dp*con2*(helpf*(-c2)+helpz*c1*nr)
      
      Edr2 = Edr2 + 2._dp*con2*(helpz*(-c2/r)+helpr*c1*nf/r)
      Edf2 = Edf2 + 2._dp*con2*c1*(helpr*(nr/r)+helpf*(2._dp*nf/r)+helpz*(nz/r))
      Edz2 = Edz2 + 2._dp*con2*(helpr*(c2/r)+helpz*nf/r)
         

  
    
    end subroutine


    subroutine egrad(n,alpha,beta,e,ga,gb)
      !!! Calculate integral free energy E and derivatives dE/da(i), dE/db(i)
      !
      ! alpha(0) - central point, alpha(i_r+(i_f-1)*nrmax) - (i_r, i_f) point
      ! NOTE THAT n vector in at origin is expressed in the direction of 1/2*df and appropriately 
      ! copied to all directions


      !THIS SUBROUTINE CALCULATES ENERGY AND DERIVATIVES OF THE ENERGY WRT GRID POINTS.
      !IT USES GAUSSIAN QUADRATURE FOR INTEGRALS AND BILINEAR INTERPOLATION IN COMPONENTS
      !OF THE N-VECTOR. 


      IMPLICIT NONE
      INTEGER :: i,ii,n,k,kk,kkk,in1,in2
      INTEGER, DIMENSION(1:4,1:2) :: ind
      REAL (KIND=dp), DIMENSION(0:nmaxt) :: alpha,beta,ga,gb
      REAL (KIND=dp), DIMENSION(0:nrmax,1:nfmax+1) :: alpha_mat,beta_mat,apsi_mat,ga_mat,gb_mat
      REAL (KIND=dp) :: bpm,apm, bc, gbc
      REAL (KIND=dp) :: apsipm,e,ap,bp,am,bm, nor
      REAL (KIND=dp) :: nrpm, nfpm, nzpm, dnr1, dnr2, dnf1, dnf2, dnz1, dnz2
     
   !   REAL (KIND=qp) :: nzpm_qp, help_qp, help2_qp, dbnz_qp  !extra precision
      
      REAL (KIND=dp) :: help, help2, epsmch
      REAL (KIND=dp), DIMENSION(1:2) :: vzpm,vrpm,vfpm,lzpm,lrpm,lfpm,wpm
     
      REAL (KIND=dp) :: da1,db1,da2,db2
      REAL (KIND=dp) :: danr_r, danr_o, danf_f, danf_o, dbnr_r, dbnr_o, dbnf_f, dbnf_o, dbnz_z, dbnz_o, danz_o
      REAL (KIND=dp) :: E0,Er,Ef,Ez, Edr1, Edf1, Edz1, Edr2, Edf2, Edz2, Ea, Eb, Eda2, Edb2
      REAL (KIND=dp),DIMENSION(1:2) :: rpm, fpm
      REAL (KIND=dp),DIMENSION(1:4,1:2,1:2) :: spm,gr1,gr2

      REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp
      REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp

      REAL (KIND=dp) :: f, ft, df, nfhelp, testhelp, help1
      nfhelp=1._dp/REAL(nfmax)
      df=2._dp*pi*nfhelp
      
      !epsmch =  EPSILON(1.0_dp) !machine precision, used to mask singularities
      
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


      
      alpha_mat=0._dp
      beta_mat=0._dp
      apsi_mat=0._dp
      ga_mat=0._dp
      gb_mat=0._dp      

      do ii=1,nfmax
         do i=1,nrmax
            alpha_mat(i,ii) = alpha(i+(ii-1)*nrmax)
            beta_mat(i,ii) =  beta(i+(ii-1)*nrmax)
            apsi_mat(i,ii) = apsi(i+(ii-1)*nrmax)
         enddo
         ft =  ii*df - df*0.5_dp 
         alpha_mat(0,ii) = alpha(0) +  (ii-1)*df  !take special status of origin into account
         beta_mat(0,ii) = beta(0)
         apsi_mat(0,ii) = apsi(0)  
      enddo
      
      !for periodic boundaries in f direction, copy first column to the nfmax+1:th column
      do i=0,nrmax         
         alpha_mat(i,nfmax+1) = alpha_mat(i,1)
         beta_mat(i,nfmax+1) = beta_mat(i,1)
         apsi_mat(i,nfmax+1) = apsi_mat(i,1)         
      enddo
      

     
     
     e=0._dp
     Ea=0._dp
     Eb=0._dp
     Eda2=0._dp     
     Edb2=0._dp
     Er=0._dp
     Ef=0._dp
     Ez=0._dp
     Edr1=0._dp
     Edr2=0._dp
     Edf1=0._dp
     Edf2=0._dp
     Edz1=0._dp
     Edz2=0._dp
   
      
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
         fpm(2)=(ii+sp)*df - df*0.5_dp
         fpm(1)=(ii+sm)*df - df*0.5_dp
         !periodic boundaries: first point at 0.5*df 
         !(phi is not needed)
         !phip=(ii+sp)*df-0.5_dp*df
         !phim=(ii+sm)*df-0.5_dp*df


      !k and kk run over interpolation grid: k in r-direction, kk in phi-direction
      do k=1,2
      do kk=1,2
        
        nrpm =0._dp 
        nfpm =0._dp
        nzpm =0._dp 
        
        dnr1 = 0._dp
        dnr2 = 0._dp
        dnf1 = 0._dp
        dnf2 = 0._dp
        dnz1 = 0._dp
        dnz2 = 0._dp
         
	apsipm=0._dp


      do kkk=1,4 !kkk runs over real grid points (corners)

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)

         !example:
         !apm =spm(1,k,kk)*alpha_mat(i,ii)+spm(2,k,kk)*alpha_mat(i+1,ii) & 
         !     +spm(3,k,kk)*alpha_mat(i,ii+1)+spm(4,k,kk)*alpha_mat(i+1,ii+1)
         apsipm = apsipm + spm(kkk,k,kk)*apsi_mat(in1,in2)
         
         !interpolate n in cylindrical component form
         nrpm = nrpm - spm(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))
         nfpm = nfpm + spm(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))
         nzpm = nzpm + spm(kkk,k,kk)*cos(beta_mat(in1,in2))
         

         dnr1 = dnr1 - gr1(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))
         dnr2 = dnr2 - gr2(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))

         dnf1 = dnf1 + gr1(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))
         dnf2 = dnf2 + gr2(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))

         dnz1 = dnz1 + gr1(kkk,k,kk)*cos(beta_mat(in1,in2))
         dnz2 = dnz2 + gr2(kkk,k,kk)*cos(beta_mat(in1,in2))



      enddo !Loop over 4 grid points (corners)
	!cut loop to include all 4 corner points at each interpolation point: en_bulk is nonlinear
      

      !check if interpolation is well defined
      nor=sqrt(nrpm**2+nzpm**2+nfpm**2)
      

      if(nor<10._dp**(-2))then
         interpolation_error=1
         WRITE(*,*) 'ERROR in interpolation: norm of n may be too small:', nor
      endif

      !normalize: divide inputs by norm     
      

      Call en_bulk(rpm(k),fpm(kk),nrpm/nor,nfpm/nor,nzpm/nor,dnr1/nor,dnr2/nor,dnf1/nor,dnf2/nor,dnz1/nor,dnz2/nor,apsipm, & 
           vzpm(k),vrpm(k),vfpm(k),lzpm(k),lrpm(k),lfpm(k), wpm(k), &
           E0,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2)
       !en_bulk call sequence:  en_bulk(r,f,nr,nf,nz,dnr1,dnr2,dnf1,dnf2,dnz1,dnz2, apsi, vz,vr,vf, & 
       !             lz,lr,lf, w, E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2)


	!write energy: not iterated over grid points
      e = e + rpm(k)*e0      

      !calculate derivatives of e wrt. grid points
      do kkk=1,4

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)   
         
         !nrpm = nrpm - spm(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))
         !nfpm = nfpm + spm(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))
         !nzpm = nzpm + spm(kkk,k,kk)*cos(beta_mat(in1,in2))

         danr_r = sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))*(nfpm**2+nzpm**2)/nor**3 
         danr_o =  -nfpm*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))/nor**3

         danf_f =  sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))*(nrpm**2+nzpm**2)/nor**3 
         danf_o =  - nrpm*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))/nor**3
         
         dbnr_r = -cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))*(nfpm**2+nzpm**2)/nor**3  
         dbnr_o = ( - nfpm*cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))  &
             + nzpm*sin(beta_mat(in1,in2)))/nor**3

         dbnf_f =  cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))*(nrpm**2+nzpm**2)/nor**3  
         dbnf_o = (nrpm*cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))  &
              + nzpm*sin(beta_mat(in1,in2)))/nor**3

         danz_o = (-nrpm*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))   &
                - nfpm*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2)))/nor**3
         dbnz_z =  - sin(beta_mat(in1,in2))*(nrpm**2+nfpm**2)/nor**3
         dbnz_o = (nrpm*cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))  &
              - nfpm*cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2)))/nor**3 
       
         
         ga_mat(in1,in2) = ga_mat(in1,in2) + Er*(danr_r+danr_o*nrpm)*spm(kkk,k,kk)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Er*(dbnr_r+dbnr_o*nrpm)*spm(kkk,k,kk)*rpm(k)
         
         ga_mat(in1,in2) = ga_mat(in1,in2) + Ef*(danf_f+danf_o*nfpm)*spm(kkk,k,kk)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Ef*(dbnf_f+dbnf_o*nfpm)*spm(kkk,k,kk)*rpm(k)

         ga_mat(in1,in2) = ga_mat(in1,in2) + Ez*(danz_o*nzpm)*spm(kkk,k,kk)*rpm(k)  !danz_z==0
         gb_mat(in1,in2) = gb_mat(in1,in2) + Ez*(dbnz_z+dbnz_o*nzpm)*spm(kkk,k,kk)*rpm(k)       
  
         
         
         !nrpm = nrpm - spm(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))
         !nfpm = nfpm + spm(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))
         !nzpm = nzpm + spm(kkk,k,kk)*cos(beta_mat(in1,in2))         

        ! dnr1 = dnr1 - gr1(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))
        ! dnr2 = dnr2 - gr2(kkk,k,kk)*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))

        ! dnf1 = dnf1 + gr1(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))
        ! dnf2 = dnf2 + gr2(kkk,k,kk)*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))

         !dnz1 = dnz1 + gr1(kkk,k,kk)*cos(beta_mat(in1,in2))
         !dnz2 = dnz2 + gr2(kkk,k,kk)*cos(beta_mat(in1,in2))


         danr_r = sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))/nor 
         danr_o =  (-nfpm*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))  & 
              -nrpm*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2) ))/nor**3

         danf_f =  sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))/nor 
         danf_o =  (-nfpm*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))  & 
              -nrpm*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2) ))/nor**3
         
         dbnr_r = -cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))/nor
         dbnr_o = ( - nfpm*cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))  &
             +nrpm*cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2)) + nzpm*sin(beta_mat(in1,in2)))/nor**3

         dbnf_f =  cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))/nor 
         dbnf_o =  ( - nfpm*cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))  &
             +nrpm*cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2)) + nzpm*sin(beta_mat(in1,in2)))/nor**3
        
         danz_o = (-nfpm*sin(beta_mat(in1,in2))*cos(alpha_mat(in1,in2))  & 
              -nrpm*sin(beta_mat(in1,in2))*sin(alpha_mat(in1,in2) ))/nor**3

         dbnz_z =  - sin(beta_mat(in1,in2))/nor
         dbnz_o = ( - nfpm*cos(beta_mat(in1,in2))*sin(alpha_mat(in1,in2))  &
             +nrpm*cos(beta_mat(in1,in2))*cos(alpha_mat(in1,in2)) + nzpm*sin(beta_mat(in1,in2)))/nor**3

         ga_mat(in1,in2) = ga_mat(in1,in2) + Edr1*(danr_r*gr1(kkk,k,kk)+danr_o*spm(kkk,k,kk)*dnr1)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edr1*(dbnr_r*gr1(kkk,k,kk)+dbnr_o*spm(kkk,k,kk)*dnr1)*rpm(k)
         
         ga_mat(in1,in2) = ga_mat(in1,in2) + Edf1*(danf_f*gr1(kkk,k,kk)+danf_o*spm(kkk,k,kk)*dnf1)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edf1*(dbnf_f*gr1(kkk,k,kk)+dbnf_o*spm(kkk,k,kk)*dnf1)*rpm(k)

         ga_mat(in1,in2) = ga_mat(in1,in2) + Edz1*(danz_o*spm(kkk,k,kk)*dnz1)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edz1*(dbnz_z*gr1(kkk,k,kk)+dbnz_o*spm(kkk,k,kk)*dnz1)*rpm(k)


         ga_mat(in1,in2) = ga_mat(in1,in2) + Edr2*(danr_r*gr2(kkk,k,kk)+danr_o*spm(kkk,k,kk)*dnr2)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edr2*(dbnr_r*gr2(kkk,k,kk)+dbnr_o*spm(kkk,k,kk)*dnr2)*rpm(k)
         
         ga_mat(in1,in2) = ga_mat(in1,in2) + Edf2*(danf_f*gr2(kkk,k,kk)+danf_o*spm(kkk,k,kk)*dnf2)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edf2*(dbnf_f*gr2(kkk,k,kk)+dbnf_o*spm(kkk,k,kk)*dnf2)*rpm(k)

         ga_mat(in1,in2) = ga_mat(in1,in2) + Edz2*(danz_o*spm(kkk,k,kk)*dnz2)*rpm(k)
         gb_mat(in1,in2) = gb_mat(in1,in2) + Edz2*(dbnz_z*gr2(kkk,k,kk)+dbnz_o*spm(kkk,k,kk)*dnz2)*rpm(k)
        

      enddo !loop over 4 grid points (corners)


      enddo !quadrature, f-direction
      enddo !quadrature, r-direction

      enddo !f-loop
      enddo !r-loop
      !set correct scale (quadrature weights and integration elements)
      e=e*0.25_dp*dx*df
      ga_mat=ga_mat*0.25_dp*dx*df
      gb_mat=gb_mat*0.25_dp*dx*df    

      
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
      e = e + E0*df*0.5_dp !no need to multiply by radius,
      !radius included in en_surf (and total radius is unity)
      ga_mat(in1,ii) = ga_mat(in1,ii) + (Ea*sm*df - Eda2)*0.5_dp
      gb_mat(in1,ii) = gb_mat(in1,ii) + (Eb*sm*df - Edb2)*0.5_dp
      ga_mat(in1,ii+1) = ga_mat(in1,ii+1) + (Ea*sp*df + Eda2)*0.5_dp
      gb_mat(in1,ii+1) = gb_mat(in1,ii+1) + (Eb*sp*df + Edb2)*0.5_dp

      call en_surf(am,bm,E0,Ea,Eb,da2,db2,Eda2,Edb2)
     
      
      !en_surf call sequence:(a,b,E,Ea,Eb,da2,db2,Eda2,Edb2)
      e = e + E0*df*0.5_dp !no need to multiply by radius, 
      !radius included in en_surf (and total radius is unity)
      ga_mat(in1,ii) = ga_mat(in1,ii) + (Ea*sp*df - Eda2)*0.5_dp
      gb_mat(in1,ii) = gb_mat(in1,ii) + (Eb*sp*df - Edb2)*0.5_dp
      ga_mat(in1,ii+1) = ga_mat(in1,ii+1) + (Ea*sm*df + Eda2)*0.5_dp
      gb_mat(in1,ii+1) = gb_mat(in1,ii+1) + (Eb*sm*df + Edb2)*0.5_dp

      enddo
      
     !CENTRAL TERM
     
    ! bc=beta_mat(0,1) !latter index irrelevant: all values are equal
    ! gbc=0._dp
    ! call en_center(bc,E0,gbc)
     
    ! e=e+E0
     


               
      ga=0._dp
      gb=0._dp
      
      !copy gradients from help slots
      do i=0,nrmax 
        ga_mat(i,1)=ga_mat(i,1)+ga_mat(i,nfmax+1)
        gb_mat(i,1)=gb_mat(i,1)+gb_mat(i,nfmax+1)                 
      enddo

      !write gradiets for testing
      
      !do i=0,nrmax
      !do ii=1,nfmax
      !   WRITE(26,*) i, ii, alpha_mat(i,ii), beta_mat(i,ii), ga_mat(i,ii), gb_mat(i,ii)         
      !enddo
      !enddo

      !reshape gradient matrices into vectors
      do ii=1,nfmax
      do i=1,nrmax
      ga(i+(ii-1)*nrmax) =  ga_mat(i,ii)
      gb(i+(ii-1)*nrmax) =  gb_mat(i,ii)     
      enddo
      ga(0) = ga(0) + ga_mat(0,ii)
      gb(0) = gb(0) + gb_mat(0,ii)
      enddo    
    ! gb(0) = gb(0) + gbc

     
    
    end subroutine

    subroutine sfun(n,x,f,g)
      !! Wrapper for egrad function for using in the TN.
      !! Calculate f and g from x values.
      !! x as array of both alpha and beta values
      !! g is array of both ga, gb
      IMPLICIT NONE
      INTEGER :: i,n !,nmax
      REAL (KIND=dp) :: f, e1, ea1, eb1, eda1, edb1, eda2, edb2, ftest
      REAL (KIND=dp), DIMENSION(n) :: x,g
      REAL (KIND=dp), DIMENSION(0:nmaxt) :: alpha,beta,ga,gb,gat,gbt, betat
      
     ! INTEGER :: mode
     ! REAL (KIND=dp), DIMENSION(n) :: gp, xp
     ! REAL (KIND=dp) :: err, fp
     ! REAL (KIND=dp), DIMENSION(0:nmaxt) :: alphap,betap,gap,gbp

      !nmax=(n-1)/2
      do i=0,nmaxt
         alpha(i)=x(i+1)
         beta(i)=x(i+nmaxt+2)
      enddo
      

      call egrad(nmaxt,alpha,beta,f,ga,gb)

      do i=0,nmaxt
         g(i+1)=ga(i)
         g(i+nmaxt+2)=gb(i) 
        ! WRITE(27,*) i, alpha(i), beta(i), ga(i), gb(i)
         
      enddo     
    
      
        
    end subroutine


END MODULE

