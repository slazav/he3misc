MODULE free

  USE modu
  USE chkder

  IMPLICIT NONE

  REAL (KIND=dp) :: s3 = sqrt(3._dp) ! tmp
  REAL (KIND=dp) :: s5 = sqrt(5._dp) ! tmp
  

  CONTAINS
    
    


    subroutine ctransf(r,f,z,zt,dzr,dzf,dzz)
      !build coordinate transformation
      
      IMPLICIT NONE
     ! INTEGER :: i
      REAL (KIND=dp) :: r,f,z,zt,dzr,dzf,dzz,h1
      
      h1=zrange
      
      zt=(r*COS(f)*TAN(stilt)/h1+1._dp)*z
      dzr = (COS(f)*TAN(stilt)/h1)*z
      dzf = (-r*SIN(f)*TAN(stilt)/h1)*z
      dzz = (r*COS(f)*TAN(stilt)/h1+1._dp)      
           
    end subroutine



    subroutine en_surf(a,b,E,Ea,Eb,da2,db2,da3,db3,Eda2,Edb2,Eda3,Edb3)
      
      !! Calculate E, dE/da, dE/db, dE/da', dE/db' at surface
      !! parameters used: de, dar,xir,lsg
      IMPLICIT NONE
     ! INTEGER :: i
      REAL (KIND=dp) :: a,b,E,Ea,Eb,Eda2,Edb2,da2,db2,da3,db3,Eda3,Edb3
      REAL (KIND=dp) :: nr,nf,nz
      REAL (KIND=dp) :: sin_a, sin_b, cos_a, cos_b, sin2b,cos2b
    
      E=0._dp
      Ea=0._dp
      Eb=0._dp
      Eda2 = 0._dp
      Edb2 = 0._dp
      Eda3=0._dp
      Edb3=0._dp

      cos_a = cos(a)
      sin_a = sin(a)
      cos_b = cos(b)
      sin_b = sin(b)
      sin2b = sin(2._dp*b)
      cos2b = cos(2._dp*b)

      nr=-sin_b*cos_a
      nf=sin_b*sin_a
      nz=cos_b

      !surface term
      E = E - 5._dp*dar*(s5*nz*nr-s3*nf)**2/16._dp
      Eb = Eb + 5._dp*dar*(s5*nz*nr-s3*nf)*(s5*cos2b*cos_a+s3*cos_b*sin_a)/8._dp
      Ea = Ea - 5._dp*dar*(s5*nz*nr-s3*nf)*(s5*nz*nf + s3*nr)/8._dp
         
      
    
      !surface gradient contribution
      E = E - 2._dp*lsg*xir**2*sin_b*(sin_b*(1-da2)+sqrt(3._dp/5._dp)*db2)/13._dp

      Eb = Eb - 2._dp*lsg*xir**2*cos_b*(2*sin_b*(1-da2)+sqrt(3._dp/5._dp)*db2)/13._dp
      
      Eda2 = Eda2 + 2._dp*lsg*xir**2*sin_b**2/13._dp
      Edb2 = Edb2 - 2._dp*lsg*xir**2*sin_b*sqrt(3._dp/5._dp)/13._dp


      !surface gradient terms, z-contr (see textcaln_3D_CORRECT.nb)

      if(1==1) then
      E = E - lsg*xir**2*da3*(2._dp*sqrt(3._dp/5._dp)*cos_a*sin_b-sin_a*sin2b)/13._dp
      E = E - 2._dp*lsg*xir**2*db3*(cos_a+sqrt(3._dp/5._dp)*cos_b*sin_a)/13._dp
      
      Ea = Ea - lsg*xir**2*(da3*(-2._dp*sqrt(3._dp/5._dp)*sin_a*sin_b-cos_a*sin2b) + &
           2._dp*db3*((-sin_a)+sqrt(3._dp/5._dp)*cos_b*cos_a))/13._dp
      Eb = Eb - lsg*xir**2*(da3*(2._dp*sqrt(3._dp/5._dp)*cos_a*cos_b-sin_a*2._dp*cos2b) + &
           2._dp*db3*(-sqrt(3._dp/5._dp)*sin_b*sin_a))/13._dp
      
      Eda3 = Eda3 - lsg*xir**2*(2._dp*sqrt(3._dp/5._dp)*cos_a*sin_b-sin_a*sin2b)/13._dp
      Edb3 = Edb3 - 2._dp*lsg*xir**2*(cos_a+sqrt(3._dp/5._dp)*cos_b*sin_a)/13._dp
      endif
      

    end subroutine

    subroutine en_freesurf(nr,nf,nz,dnr1,dnr2,dnf1,dnf2,dnz1,dnz2,dnr3,dnf3,dnz3, & 
         E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2,Edr3,Edf3,Edz3)      

      REAL (KIND=dp) :: r,f,E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2     
      REAL (KIND=dp) :: nz,nr,nf, rzz,rzr,rzf, dnr1,dnr2,dnf1,dnf2,dnz1,dnz2
      REAL (KIND=dp) :: dnr3, dnf3, dnz3, Edr3, Edf3, Edz3
      REAL (KIND=dp) ::c,s


      c=-0.25_dp !\cos\theta 
      s=SQRT(15.)/4.0_dp !\sin\theta 
      rzr=(1-c)*nz*nr-s*nf ! H*Rij
      rzf=(1-c)*nz*nf+s*nr
      rzz=c+(1-c)*nz**2


    end subroutine en_freesurf


    subroutine en_bulk(r,f,nr,nf,nz,dnr1,dnr2,dnf1,dnf2,dnz1,dnz2, apsi, & 
         vz,vr,vf, lz,lr,lf, w, dnr3,dnf3,dnz3, & 
         E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2, Edr3,Edf3,Edz3)      
      
      !! parameters used:
      !!   chia*(nub/nu0)^2, for non-zero apsi
      !!   lo for non-zero rotation
      !!   vd for non-zero flow
      !!   de and xir
      REAL (KIND=dp) :: r,f,E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2
      REAL (KIND=dp) :: apsi, vz,vr,vf, lz,lr,lf, w
      REAL (KIND=dp) :: nz,nr,nf, rzz,rzr,rzf, dnr1,dnr2,dnf1,dnf2,dnz1,dnz2
      REAL (KIND=dp) :: dnr3, dnf3, dnz3, Edr3, Edf3, Edz3
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

      Edr3 = 0._dp
      Edf3 = 0._dp
      Edz3 = 0._dp


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
      !See: .../spinwaves/.../Derivations/textcalcn_2D_CORRECT.nb (see wikipedia for diff operators in cylindrical coord)
      !Here we use the tensor expression directly (transformed into vector component form), see handwritten notes (Samuli)


      !see textcalcn_2D_CORRECT.nb
      !con1 = 2._dp*(4._dp+de)*xir**2/13._dp      
      con1 = 8._dp*xir**2/13._dp    !see Kopus original 1D (with slazav modifications for readibility)
      !WRITE(*,*)con1

      ! con1 = 8._dp*(9._dp+2._dp*de)*xir**2/65._dp
      !! con1= 0.5_dp*xir**2
      


      !(\nabla_i n_j)^2, cyl coord, see textcal_3D_nablai_nj.nb
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

      !(\nabla_i n_j)^2 (z-derivatives), cyl coord (see textcaln_3D_correct_nablai_nj.nb)
       
      E = E + con1*(dnr3**2+dnf3**2+dnz3**2)
      
      Edr3 = Edr3 + 2._dp*con1*dnr3
      Edf3 = Edf3 + 2._dp*con1*dnf3
      Edz3 = Edz3 + 2._dp*con1*dnz3
      
      
      
      
      
      
      con2 = 8._dp*(2._dp+de)*xir**2/65._dp  !see textcalcn_2D_CORRECT.nb and Kopus original 1D (with perhaps slazav modifications for readibility)
      !!con2 = xir**2/5._dp
      !WRITE(*,*)con2
      ! {5/4*(n /div n +(n /nabla) n)-sqrt(15)/4 /rot(n)}^2 , see "material derivative" for (n /nabla)

      helpr = c1*(nr*(dnr1+dnf2/r+nr/r+dnz3) + nr*dnr1+nf*dnr2/r-nf**2/r + nz*dnr3) + c2*(dnz2/r -dnf3)
      helpf = c1*(nf*(dnr1+dnf2/r+nr/r+dnz3) + nr*dnf1+nf*dnf2/r+nr*nf/r +nz*dnf3 ) + c2*(-dnz1+dnr3) 
      helpz = c1*(nz*(dnr1+dnf2/r+nr/r+dnz3) + nr*dnz1+nf*dnz2/r+ nz*dnz3 ) + c2*(nf+dnf1*r-dnr2)/r 

      E = E +con2*(helpr**2+helpf**2+helpz**2)
         

      !Er = Er + 2._dp*con2*c1*(helpr*(2._dp*dnr1+dnf2/r+2._dp*nr/r)+helpf*(2._dp*nf/r+dnf1)+helpz*(nz/r+dnz1))
      !Ef = Ef + 2._dp*con2*(helpr*c1*(dnr2/r-2._dp*nf/r)+helpf*c1*(dnr1+2._dp*dnf2/r+2._dp*nr/r) +helpz*(c1*dnz2/r+c2/r))
      !Ez = Ez + 2._dp*con2*helpz*c1*(dnr1+dnf2/r+nr/r)

      
      Er = Er + 2._dp*con2*c1*(helpr*(2._dp*dnr1+dnf2/r+2._dp*nr/r+dnz3)+helpf*(2._dp*nf/r+dnf1)+helpz*(nz/r+dnz1))
      Ef = Ef + 2._dp*con2*(helpr*c1*(dnr2/r-2._dp*nf/r+dnz3)+helpf*c1*(dnr1+2._dp*dnf2/r+2._dp*nr/r) +helpz*(c1*dnz2/r+c2/r))
      Ez = Ez + 2._dp*con2*c1*(helpz*(dnr1+dnf2/r+nr/r+2._dp*dnz3)+helpr*dnr3+helpf*dnf3)

      Edr1 = Edr1 + 2._dp*con2*c1*(helpr*2._dp*nr+helpf*nf+helpz*nz)
      Edf1 = Edf1 + 2._dp*con2*(helpz*c2+c1*helpf*nr)
      Edz1 = Edz1 + 2._dp*con2*(helpf*(-c2)+helpz*c1*nr)
      
      Edr2 = Edr2 + 2._dp*con2*(helpz*(-c2/r)+helpr*c1*nf/r)
      Edf2 = Edf2 + 2._dp*con2*c1*(helpr*(nr/r)+helpf*(2._dp*nf/r)+helpz*(nz/r))
      Edz2 = Edz2 + 2._dp*con2*(helpr*(c2/r)+helpz*nf/r)
         
      Edr3 = Edr3 + 2._dp*con2*(c1*helpr*nz+c2*helpf)
      Edf3 = Edf3 + 2._dp*con2*(c1*helpf*nz-c2*helpr)
      Edz3 = Edz3 + 2._dp*c1*con2*(helpr*nr+helpf*nf+2._dp*helpz*nz)
  
    
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
      INTEGER :: i,ii,n,k,kk,kkk,in1,in2,in3, iii, k3
      INTEGER, DIMENSION(1:8,1:3) :: ind
      INTEGER, DIMENSION(1:4,1:2) :: sind
      REAL (KIND=dp), DIMENSION(0:nmaxt-1) :: alpha,beta,ga,gb
      REAL (KIND=dp), DIMENSION(0:nrmax,1:nfmax+1,1:nzmax) :: alpha_mat,beta_mat,apsi_mat,ga_mat,gb_mat
      REAL (KIND=dp) :: bpm,apm, bc, gbc
      REAL (KIND=dp) :: apsipm,e,ap,bp,am,bm, nor
      REAL (KIND=dp) :: nrpm, nfpm, nzpm, dnr1, dnr2, dnf1, dnf2, dnz1, dnz2
     
   !   REAL (KIND=qp) :: nzpm_qp, help_qp, help2_qp, dbnz_qp  !extra precision
      
      REAL (KIND=dp) :: help, help2, epsmch
      REAL (KIND=dp), DIMENSION(1:2) :: vzpm,vrpm,vfpm,lzpm,lrpm,lfpm,wpm
     
      REAL (KIND=dp) :: da1,db1,da2,db2,da3,db3,zt,dzr,dzf,dzz
      REAL (KIND=dp) :: danr_r, danr_o, danf_f, danf_o, dbnr_r, dbnr_o, dbnf_f, dbnf_o, dbnz_z, dbnz_o, danz_o
      REAL (KIND=dp) :: E0,Er,Ef,Ez, Edr1, Edf1, Edz1, Edr2, Edf2, Edz2, Ea, Eb, Eda2, Edb2, Edb3, Eda3
      REAL (KIND=dp) :: dnr3, dnf3, dnz3, Edr3, Edf3, Edz3

      REAL (KIND=dp),DIMENSION(1:2) :: rpm, fpm, zpm
      REAL (KIND=dp),DIMENSION(1:8,1:2,1:2,1:2) :: spm,gr1,gr2,gr3,gr1t,gr2t,gr3t
      REAL (KIND=dp),DIMENSION(1:4,1:2,1:2) :: sspm,sgr1,sgr2, sgr3

      REAL (KIND=dp) :: sp = (3._dp + sqrt(3._dp))/6._dp
      REAL (KIND=dp) :: sm = (3._dp - sqrt(3._dp))/6._dp

      REAL (KIND=dp) :: f, ft, df, dz, nfhelp, testhelp, help1
      nfhelp=1._dp/REAL(nfmax)
      df=2._dp*pi*nfhelp
      
      dz=zrange/r/(nzmax-1._dp)  !normalize to unit radius
      

      !epsmch =  EPSILON(1.0_dp) !machine precision, used to mask singularities
      
      !periodic boundary conditions->n_intervals=n_points
      !construct help arrays

      !help indices for r, f and z directions in interpolation:
      !to take steps over corners (real grid points)
      ind(:,1)= (/ 0,1,0,1 , 0,1,0,1 /) !r-direction
      ind(:,2)= (/ 0,0,1,1 , 0,0,1,1 /) !f-direction
      ind(:,3)= (/ 0,0,0,0 , 1,1,1,1 /) !z-direction


      !spm(:,1,1) = (/ sp*sp,sp*sm,sp*sm,sm*sm /) !to first interp. point from 4 corner points
      !spm(:,2,1) = (/ sp*sm,sp*sp,sm*sm,sp*sm /) !second interp. point
      !spm(:,1,2) = (/ sp*sm,sm*sm,sp*sp,sp*sm /) !third  interp. point
      !spm(:,2,2) = (/ sm*sm,sm*sp,sm*sp,sp*sp /) !fourth interp. point

      !weights for trilinear interpolation to gaussian quadrature points (8)
      spm(1:4,1,1,1) = (/ sp*sp,sp*sm,sp*sm,sm*sm /) !to first interp. point from 4 first corner points
      spm(1:4,2,1,1) = (/ sp*sm,sp*sp,sm*sm,sp*sm /) !second interp. point
      spm(1:4,1,2,1) = (/ sp*sm,sm*sm,sp*sp,sp*sm /) !third  interp. point
      spm(1:4,2,2,1) = (/ sm*sm,sm*sp,sm*sp,sp*sp /) !fourth interp. point

      spm(5:8,1:2,1:2,1) = spm(1:4,1:2,1:2,1)*sm 
      spm(1:4,1:2,1:2,2) = spm(1:4,1:2,1:2,1)*sm
      spm(5:8,1:2,1:2,2) = spm(1:4,1:2,1:2,1)*sp
      spm(1:4,1:2,1:2,1) = spm(1:4,1:2,1:2,1)*sp
      
      gr1t(1:4,1,1,1) = (/ -sp,sp,-sm,sm /) !first interp. point from 4 corner points
      gr1t(1:4,2,1,1) = (/ -sp,sp,-sm,sm /) !second interp. point
      gr1t(1:4,1,2,1) = (/ -sm,sm,-sp,sp /) !third interp. point
      gr1t(1:4,2,2,1) = (/ -sm,sm,-sp,sp /) !fourth interp. point
      
      gr1t(5:8,1:2,1:2,1) = gr1t(1:4,1:2,1:2,1)*sm
      gr1t(1:4,1:2,1:2,2) = gr1t(1:4,1:2,1:2,1)*sm
      gr1t(5:8,1:2,1:2,2) = gr1t(1:4,1:2,1:2,1)*sp
      gr1t(1:4,1:2,1:2,1) = gr1t(1:4,1:2,1:2,1)*sp

      gr1t=gr1t/dx
     


      gr2t(1:4,1,1,1) = (/ -sp,-sm,sp,sm /) !first interp. point
      gr2t(1:4,2,1,1) = (/ -sm,-sp,sm,sp /) !second interp. point
      gr2t(1:4,1,2,1) = (/ -sp,-sm,sp,sm /) !third interp. point
      gr2t(1:4,2,2,1) = (/ -sm,-sp,sm,sp /) !fourth interp. point

      gr2t(5:8,1:2,1:2,1) = gr2t(1:4,1:2,1:2,1)*sm
      gr2t(1:4,1:2,1:2,2) = gr2t(1:4,1:2,1:2,1)*sm
      gr2t(5:8,1:2,1:2,2) = gr2t(1:4,1:2,1:2,1)*sp
      gr2t(1:4,1:2,1:2,1) = gr2t(1:4,1:2,1:2,1)*sp
      gr2t=gr2t/df


      gr3t(1:4,1,1,1) = (/ sp*sp,sp*sm,sp*sm,sm*sm /) !to first interp. point from 4 corner points
      gr3t(1:4,2,1,1) = (/ sp*sm,sp*sp,sm*sm,sp*sm /) !second interp. point
      gr3t(1:4,1,2,1) = (/ sp*sm,sm*sm,sp*sp,sp*sm /) !third  interp. point
      gr3t(1:4,2,2,1) = (/ sm*sm,sm*sp,sm*sp,sp*sp /) !fourth interp. point

      gr3t(5:8,1:2,1:2,1) = gr3t(1:4,1:2,1:2,1)*(-1._dp)
      gr3t(1:4,1:2,1:2,2) = gr3t(1:4,1:2,1:2,1)*(-1._dp)
      gr3t(5:8,1:2,1:2,2) = gr3t(1:4,1:2,1:2,1)
      gr3t(1:4,1:2,1:2,1) = gr3t(1:4,1:2,1:2,1)

      gr3t=gr3t/dz
      
          
    

      alpha_mat=0._dp
      beta_mat=0._dp
      apsi_mat=0._dp
      ga_mat=0._dp
      gb_mat=0._dp      

      
      
   do iii=1,nzmax
      do ii=1,nfmax
         do i=1,nrmax
            alpha_mat(i,ii,iii) = alpha(i+(ii-1)*nrmax+(iii-1)*nmaxp)
            beta_mat(i,ii,iii) =  beta(i+(ii-1)*nrmax+(iii-1)*nmaxp)
            apsi_mat(i,ii,iii) = apsi(i+(ii-1)*nrmax+(iii-1)*nmaxp)
         enddo
         ft =  ii*df - df*0.5_dp 
         alpha_mat(0,ii,iii) = alpha(0+(iii-1)*nmaxp) +  (ii-0.5_dp)*df  !take special status of origin into account
         beta_mat(0,ii,iii) = beta(0+(iii-1)*nmaxp)
         apsi_mat(0,ii,iii) = apsi(0+(iii-1)*nmaxp)  
      enddo
      

      !for periodic boundaries in f direction, copy first column to the nfmax+1:th column
      do i=0,nrmax         
         alpha_mat(i,nfmax+1,iii) = alpha_mat(i,1,iii)
         beta_mat(i,nfmax+1,iii) = beta_mat(i,1,iii)
         apsi_mat(i,nfmax+1,iii) = apsi_mat(i,1,iii)         
      enddo
      
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
     Edr3=0._dp
     Edf3=0._dp
     Edz3=0._dp
   
      
      do i=0,nrmax-1
         !f,z-independent iterpolation
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

      !f-dependent,z-independent interpolation
      do ii=1,nfmax
         fpm(2)=(ii+sp)*df - df*0.5_dp
         fpm(1)=(ii+sm)*df - df*0.5_dp
         !periodic boundaries: first point at 0.5*df 
         !(phi is not needed)
         !phip=(ii+sp)*df-0.5_dp*df
         !phim=(ii+sm)*df-0.5_dp*df
      


     ! z-dependent interpolation
      do iii=1,nzmax-1
         zpm(1)=(iii+sp)*dz
         zpm(2)=(iii+sm)*dz

     
      


      !k and kk run over interpolation grid: k in r-direction, kk in phi-direction
      do k=1,2
      do kk=1,2
      do k3=1,2

        !inverse coordinate transformation call sequence: ctransf(r,f,z,zt,dzr,dzf,dzz)
        CALL ctransf(rpm(k),fpm(kk),zpm(k3),zt,dzr,dzf,dzz)
        gr1=gr1t+gr3t*dzr
        gr2=gr2t+gr3t*dzf
        gr3=gr3t*dzz


        nrpm =0._dp 
        nfpm =0._dp
        nzpm =0._dp 
        
        dnr1 = 0._dp
        dnr2 = 0._dp
        dnf1 = 0._dp
        dnf2 = 0._dp
        dnz1 = 0._dp
        dnz2 = 0._dp
        dnr3 = 0._dp
        dnf3 = 0._dp
        dnz3 = 0._dp
         
	apsipm=0._dp


      do kkk=1,8 !kkk runs over real grid points (corners)

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)
         in3=iii+ind(kkk,3)

 

         !example:
         !apm =spm(1,k,kk)*alpha_mat(i,ii)+spm(2,k,kk)*alpha_mat(i+1,ii) & 
         !     +spm(3,k,kk)*alpha_mat(i,ii+1)+spm(4,k,kk)*alpha_mat(i+1,ii+1)
        
         apsipm = apsipm + spm(kkk,k,kk,k3)*apsi_mat(in1,in2,in3)
         
         !interpolate n in cylindrical component form
         nrpm = nrpm - spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         nfpm = nfpm + spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         nzpm = nzpm + spm(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))
         

         dnr1 = dnr1 - gr1(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         dnr2 = dnr2 - gr2(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         dnr3 = dnr3 - gr3(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))

         dnf1 = dnf1 + gr1(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         dnf2 = dnf2 + gr2(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         dnf3 = dnf3 + gr3(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))

         dnz1 = dnz1 + gr1(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))
         dnz2 = dnz2 + gr2(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))
         dnz3 = dnz3 + gr3(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))


        

      enddo !Loop over 8 grid points (corners)
	!cut loop to include all 8 corner points at each interpolation point: en_bulk is nonlinear
      

      !check if interpolation is well defined
      nor=sqrt(nrpm**2+nzpm**2+nfpm**2)
      

      if(nor<10._dp**(-2))then
         interpolation_error=1
         WRITE(*,*) 'ERROR in interpolation: norm of n may be too small:', nor
      endif

      !normalize: divide inputs by norm     
       

      Call en_bulk(rpm(k),fpm(kk),nrpm/nor,nfpm/nor,nzpm/nor,dnr1/nor,dnr2/nor,dnf1/nor,dnf2/nor,dnz1/nor,dnz2/nor,apsipm, & 
           vzpm(k),vrpm(k),vfpm(k),lzpm(k),lrpm(k),lfpm(k), wpm(k), dnr3/nor,dnf3/nor,dnz3/nor, &
           E0,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2, Edr3,Edf3,Edz3)
       !en_bulk call sequence:  subroutine en_bulk(r,f,nr,nf,nz,dnr1,dnr2,dnf1,dnf2,dnz1,dnz2, apsi, & 
      !                                      vz,vr,vf, lz,lr,lf, w, dnr3,dnf3,dnz3, & 
       !                                E,Er,Ef,Ez,Edr1,Edr2,Edf1,Edf2,Edz1,Edz2, Edr3,Edf3,Edz3)     

	!write energy: not iterated over grid points
      e = e + rpm(k)*e0*ABS(dzz)


     

      !calculate derivatives of e wrt. grid points
      do kkk=1,8

         in1=i+ind(kkk,1)
         in2=ii+ind(kkk,2)   
         in3=iii+ind(kkk,3)

         !nrpm = nrpm - spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         !nfpm = nfpm + spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         !nzpm = nzpm + spm(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))

         danr_r = sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))*(nfpm**2+nzpm**2)/nor**3 
         danr_o =  -nfpm*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))/nor**3

         danf_f =  sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))*(nrpm**2+nzpm**2)/nor**3 
         danf_o =  - nrpm*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))/nor**3
         
         dbnr_r = -cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))*(nfpm**2+nzpm**2)/nor**3  
         dbnr_o = ( - nfpm*cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))  &
             + nzpm*sin(beta_mat(in1,in2,in3)))/nor**3

         dbnf_f =  cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))*(nrpm**2+nzpm**2)/nor**3  
         dbnf_o = (nrpm*cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))  &
              + nzpm*sin(beta_mat(in1,in2,in3)))/nor**3

         danz_o = (-nrpm*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))   &
                - nfpm*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3)))/nor**3
         dbnz_z =  - sin(beta_mat(in1,in2,in3))*(nrpm**2+nfpm**2)/nor**3
         dbnz_o = (nrpm*cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))  &
              - nfpm*cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3)))/nor**3 
       
         
         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Er*(danr_r+danr_o*nrpm)*spm(kkk,k,kk,k3)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Er*(dbnr_r+dbnr_o*nrpm)*spm(kkk,k,kk,k3)*rpm(k)
         
         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Ef*(danf_f+danf_o*nfpm)*spm(kkk,k,kk,k3)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Ef*(dbnf_f+dbnf_o*nfpm)*spm(kkk,k,kk,k3)*rpm(k)

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Ez*(danz_o*nzpm)*spm(kkk,k,kk,k3)*rpm(k)  !danz_z==0
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Ez*(dbnz_z+dbnz_o*nzpm)*spm(kkk,k,kk,k3)*rpm(k)       
  
         
         
   
       
         !nrpm = nrpm - spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         !nfpm = nfpm + spm(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         !nzpm = nzpm + spm(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))


         !dnr1 = dnr1 - gr1(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         !dnr2 = dnr2 - gr2(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))
         !dnr3 = dnr3 - gr3(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))

         !dnf1 = dnf1 + gr1(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         !dnf2 = dnf2 + gr2(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))
         !dnf3 = dnf3 + gr3(kkk,k,kk,k3)*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))

         !dnz1 = dnz1 + gr1(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))
         !dnz2 = dnz2 + gr2(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))
         !dnz3 = dnz3 + gr3(kkk,k,kk,k3)*cos(beta_mat(in1,in2,in3))


         danr_r = sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))/nor 
         danr_o =  (-nfpm*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))  & 
              -nrpm*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3) ))/nor**3

         danf_f =  sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))/nor 
         danf_o =  (-nfpm*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))  & 
              -nrpm*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3) ))/nor**3
         
         dbnr_r = -cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))/nor
         dbnr_o = ( - nfpm*cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))  &
             +nrpm*cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3)) + nzpm*sin(beta_mat(in1,in2,in3)))/nor**3

         dbnf_f =  cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))/nor 
         dbnf_o =  ( - nfpm*cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))  &
             +nrpm*cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3)) + nzpm*sin(beta_mat(in1,in2,in3)))/nor**3
        
         danz_o = (-nfpm*sin(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3))  & 
              -nrpm*sin(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3) ))/nor**3

         dbnz_z =  - sin(beta_mat(in1,in2,in3))/nor
         dbnz_o = ( - nfpm*cos(beta_mat(in1,in2,in3))*sin(alpha_mat(in1,in2,in3))  &
             +nrpm*cos(beta_mat(in1,in2,in3))*cos(alpha_mat(in1,in2,in3)) + nzpm*sin(beta_mat(in1,in2,in3)))/nor**3

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edr1*(danr_r*gr1(kkk,k,kk,k3)+danr_o*spm(kkk,k,kk,k3)*dnr1)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edr1*(dbnr_r*gr1(kkk,k,kk,k3)+dbnr_o*spm(kkk,k,kk,k3)*dnr1)*rpm(k)
         
         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edf1*(danf_f*gr1(kkk,k,kk,k3)+danf_o*spm(kkk,k,kk,k3)*dnf1)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edf1*(dbnf_f*gr1(kkk,k,kk,k3)+dbnf_o*spm(kkk,k,kk,k3)*dnf1)*rpm(k)

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edz1*(danz_o*spm(kkk,k,kk,k3)*dnz1)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edz1*(dbnz_z*gr1(kkk,k,kk,k3)+dbnz_o*spm(kkk,k,kk,k3)*dnz1)*rpm(k)


         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edr2*(danr_r*gr2(kkk,k,kk,k3)+danr_o*spm(kkk,k,kk,k3)*dnr2)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edr2*(dbnr_r*gr2(kkk,k,kk,k3)+dbnr_o*spm(kkk,k,kk,k3)*dnr2)*rpm(k)
         
         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edf2*(danf_f*gr2(kkk,k,kk,k3)+danf_o*spm(kkk,k,kk,k3)*dnf2)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edf2*(dbnf_f*gr2(kkk,k,kk,k3)+dbnf_o*spm(kkk,k,kk,k3)*dnf2)*rpm(k)

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edz2*(danz_o*spm(kkk,k,kk,k3)*dnz2)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edz2*(dbnz_z*gr2(kkk,k,kk,k3)+dbnz_o*spm(kkk,k,kk,k3)*dnz2)*rpm(k)
        

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edr3*(danr_r*gr3(kkk,k,kk,k3)+danr_o*spm(kkk,k,kk,k3)*dnr3)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edr3*(dbnr_r*gr3(kkk,k,kk,k3)+dbnr_o*spm(kkk,k,kk,k3)*dnr3)*rpm(k)
         
         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edf3*(danf_f*gr3(kkk,k,kk,k3)+danf_o*spm(kkk,k,kk,k3)*dnf3)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edf3*(dbnf_f*gr3(kkk,k,kk,k3)+dbnf_o*spm(kkk,k,kk,k3)*dnf3)*rpm(k)

         ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + Edz3*(danz_o*spm(kkk,k,kk,k3)*dnz3)*rpm(k)
         gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + Edz3*(dbnz_z*gr3(kkk,k,kk,k3)+dbnz_o*spm(kkk,k,kk,k3)*dnz3)*rpm(k)
        

      enddo !loop over 8 grid points (corners)

      enddo !quadrature, z-direction
      enddo !quadrature, f-direction
      enddo !quadrature, r-direction
      
      enddo !r-loop
      enddo !f-loop
      enddo !z-loop
     
      
      !set correct scale (quadrature weights and integration elements)
      e=e*0.125_dp*dx*df*dz
      ga_mat=ga_mat*0.125_dp*dx*df*dz
      gb_mat=gb_mat*0.125_dp*dx*df*dz    
  
      
      
 !help indices for r and f directions in interpolation:
      !to take steps over corners (real grid points)
      sind(:,1)= (/ 0,1,0,1 /) !f-direction
      sind(:,2)= (/ 0,0,1,1 /) !z-direction
      



      !weights for bilinear interpolation to gaussian quadrature points (4)
      sspm(:,1,1) = (/ sp*sp,sp*sm,sp*sm,sm*sm /) !to first interp. point from 4 corner points
      sspm(:,2,1) = (/ sp*sm,sp*sp,sm*sm,sp*sm /) !second interp. point
      sspm(:,1,2) = (/ sp*sm,sm*sm,sp*sp,sp*sm /) !third  interp. point
      sspm(:,2,2) = (/ sm*sm,sm*sp,sm*sp,sp*sp /) !fourth interp. point
      
      !weights for gradients of the interpolation (1 is f direction,2 is z direction)
      sgr2(:,1,1) = (/ -sp,sp,-sm,sm /) !first interp. point from 4 corner points
      sgr2(:,2,1) = (/ -sp,sp,-sm,sm /) !second interp. point
      sgr2(:,1,2) = (/ -sm,sm,-sp,sp /) !third interp. point
      sgr2(:,2,2) = (/ -sm,sm,-sp,sp /) !fourth interp. point
      sgr2=sgr2/df
     
      sgr3(:,1,1) = (/ -sp,-sm,sp,sm /) !first interp. point
      sgr3(:,2,1) = (/ -sm,-sp,sm,sp /) !second interp. point
      sgr3(:,1,2) = (/ -sp,-sm,sp,sm /) !third interp. point
      sgr3(:,2,2) = (/ -sm,-sp,sm,sp /) !fourth interp. point
      sgr3=sgr3/dz

      
      


      ! surface terms: NOTE THAT RADIUS IS INCLUDED IN en_surf, no need to multiply by r
      !furthermore, Rmax=1 in this code
         


     

      if(1==1) then

      in1=nrmax !boundary

      

      do iii=1,nzmax-1
      do ii=1,nfmax
      do k=1,2
      do kk=1,2

      apm=0._dp
      bpm=0._dp
      da2=0._dp
      db2=0._dp
      da3=0._dp
      db3=0._dp

      do kkk=1,4 !kkk runs over real grid points (corners)
         
        
         in2=ii+sind(kkk,1)
         in3=iii+sind(kkk,2)

        apm = apm + sspm(kkk,k,kk)*alpha_mat(in1,in2,in3)	
        
        bpm = bpm + sspm(kkk,k,kk)*beta_mat(in1,in2,in3)
        
	da2 = da2 + sgr2(kkk,k,kk)*alpha_mat(in1,in2,in3)	
	db2 = db2 + sgr2(kkk,k,kk)*beta_mat(in1,in2,in3)
      
        da3 = da3 + sgr3(kkk,k,kk)*alpha_mat(in1,in2,in3)	
	db3 = db3 + sgr3(kkk,k,kk)*beta_mat(in1,in2,in3)

     enddo !cut loop before calculating the energy

     

      call en_surf(apm,bpm,E0,Ea,Eb,da2,db2,da3,db3,Eda2,Edb2,Eda3,Edb3)      
      !en_surf call sequence:(a,b,E,Ea,Eb,da2,db2,da3,db3,Eda2,Edb2,Eda3,Edb3)
     


      e = e + E0*df*dz*0.25_dp !no need to multiply by radius
      !radius included in en_surf (and total radius is unity)
     do kkk=1,4
         in2=ii+sind(kkk,1)
         in3=iii+sind(kkk,2)

      ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + (Ea*sspm(kkk,k,kk) + Eda2*sgr2(kkk,k,kk))*0.25_dp*dz*df
      gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + (Eb*sspm(kkk,k,kk) + Edb2*sgr2(kkk,k,kk))*0.25_dp*dz*df
     
      ga_mat(in1,in2,in3) = ga_mat(in1,in2,in3) + (Eda3*sgr3(kkk,k,kk))*0.25_dp*dz*df
      gb_mat(in1,in2,in3) = gb_mat(in1,in2,in3) + (Edb3*sgr3(kkk,k,kk))*0.25_dp*dz*df

     enddo
     enddo
     enddo
     enddo
     enddo

     endif
     


               
      ga=0._dp
      gb=0._dp
      
      !copy gradients from help slots
    do iii=1,nzmax
      do i=0,nrmax 
        ga_mat(i,1,iii)=ga_mat(i,1,iii)+ga_mat(i,nfmax+1,iii)
        gb_mat(i,1,iii)=gb_mat(i,1,iii)+gb_mat(i,nfmax+1,iii)                 
      enddo
   enddo
      !write gradiets for testing
      do iii=1,nzmax
      do i=0,nrmax
      do ii=1,nfmax
         WRITE(26,*) i, ii, iii, alpha_mat(i,ii,iii), beta_mat(i,ii,iii), ga_mat(i,ii,iii), gb_mat(i,ii,iii)         
      enddo
      enddo
      enddo

      !reshape gradient matrices into vectors
    
      do iii=1,nzmax
      do ii=1,nfmax
      do i=1,nrmax
      ga(i+(ii-1)*nrmax+(iii-1)*nmaxp) =  ga_mat(i,ii,iii)
      gb(i+(ii-1)*nrmax+(iii-1)*nmaxp) =  gb_mat(i,ii,iii)     
      enddo
      ga(0+(iii-1)*nmaxp) = ga(0+(iii-1)*nmaxp) + ga_mat(0,ii,iii)
      gb(0+(iii-1)*nmaxp) = gb(0+(iii-1)*nmaxp) + gb_mat(0,ii,iii)
      enddo    
      enddo    

     
    end subroutine

    subroutine sfun(n,x,f,g)
      !! Wrapper for egrad function for using in the TN.
      !! Calculate f and g from x values.
      !! x as array of both alpha and beta values
      !! g is array of both ga, gb
      IMPLICIT NONE
      INTEGER :: i,n !,nmax
      REAL (KIND=dp) :: f, e1, ea1, eb1, eda1, edb1, eda2, edb2, ftest
      REAL (KIND=dp), DIMENSION(1:n) :: x,g
      REAL (KIND=dp), DIMENSION(0:nmaxt-1) :: alpha,beta,ga,gb,gat,gbt, betat
      
     ! INTEGER :: mode
     ! REAL (KIND=dp), DIMENSION(n) :: gp, xp
     ! REAL (KIND=dp) :: err, fp
     ! REAL (KIND=dp), DIMENSION(0:nmaxp) :: alphap,betap,gap,gbp

      

      do i=0,nmaxt-1
         alpha(i)=x(i+1)
         beta(i)=x(i+nmaxt+1)
      enddo
     

      call egrad(nmaxt,alpha,beta,f,ga,gb)

     

      do i=0,nmaxt-1
         g(i+1)=ga(i)
         g(i+nmaxt+1)=gb(i) 
         WRITE(27,*) i, alpha(i), beta(i), ga(i), gb(i)         
      enddo     
    
      
      
        
    end subroutine


END MODULE

