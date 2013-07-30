
      double precision function energy(n,m,a)
      implicit none 
      integer i,j,m,n,k,l,nmax,ilimit
      double precision pii,r,dr,df,rp
      parameter (nmax=200,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision omega0,de,d0,apu,eps,fsg,co
      double precision fh,fv,fs,fg,l0,fl,ov0,y(nmax)
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      fh=0d0
      fv=0d0
      fg=0d0
      fs=0d0
      fl=0d0
      fsg=0d0
      eps=1d-8
      rmax=dble(n-1)*dr
      ilimit=dint(dble(n)*dsqrt(ov0/omega0)+0.5d0)
      if (omega0.eq.0d0) ilimit=0
      do 200 i=2,n
         r=dr*dble(i-1)
         do 100 j=2,m-1
            if (i.gt.ilimit) then
               v=omega0*r-ov0*rmax*rmax/r
            else
               v=0d0
               fl=fl+l0*r*y(i)*a(3,3,i,j)**2
            end if
            fh=fh-r*(4d0*a(3,3,i,j)+1d0)*y(i)
            fv=fv-2d0*r*y(i)*(a(3,2,i,j)*v)**2
            fg=fg+(1d0/r)*y(i)*(4d0+de)*
     $           ((a(1,2,i,j)+a(2,1,i,j))**2+
     $           (a(1,1,i,j)-a(2,2,i,j))**2)
            fg=fg+(1d0/r)*y(i)*((3d0+de)*
     $           a(3,1,i,j)**2+a(1,3,i,j)**2+
     $           a(3,2,i,j)**2+a(2,3,i,j)**2)
 100     continue
 200  continue
      do 300 j=2,m-1
         fs=fs-5d0*(rmax/dr)*d0*a(3,1,n,j)**2
 300  continue
      do 800 i=1,n-1
         do 700 j=2,m-1
            rp=(dble(i)-0.5d0)*dr
            do 600 k=1,3
               do 500 l=1,3
                  if (l.eq.1) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  fg=fg+apu*rp*((a(k,l,i+1,j)-a(k,l,i,j))/dr)**2
 500           continue
 600        continue
 700     continue
 800  continue
      do 1000 i=2,n-1
         do 900 j=2,m-1
            fg=fg+(2d0+de)*(a(1,1,i+1,j)+a(1,1,i,j)
     $           -a(2,2,i+1,j)-a(2,2,i,j))*
     $           (a(1,1,i+1,j)-a(1,1,i,j))/dr
            fg=fg+(2d0+de)*(a(1,2,i+1,j)+a(1,2,i,j)
     $           +a(2,1,i+1,j)+a(2,1,i,j))*
     $           (a(2,1,i+1,j)-a(2,1,i,j))/dr
            fg=fg+(2d0+de)*(a(3,1,i+1,j)+a(3,1,i,j))
     $           *(a(3,1,i+1,j)-a(3,1,i,j))/dr
 900     continue
 1000 continue
      do 1700 i=2,n
c         do 1600 j=1,m-2
         do 1600 j=2,m-2
            r=dble(i-1)*dr
            do 1500 k=1,3
               do 1400 l=1,3
                  if (l.eq.3) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  fg=fg+(apu/r)*((a(k,l,i,j+1)-a(k,l,i,j))/df)**2
 1400          continue
 1500       continue
            fg=fg+(3d0+de)*(a(1,1,i,j+1)+a(1,1,i,j)
     $           -a(2,2,i,j+1)-a(2,2,i,j))*
     $           (a(1,2,i,j+1)-a(1,2,i,j))/(df*r)
            fg=fg+(3d0+de)*(a(1,2,i,j+1)+a(1,2,i,j)
     $           +a(2,1,i,j+1)+a(2,1,i,j))*
     $           (a(2,2,i,j+1)-a(2,2,i,j))/(df*r)
            fg=fg+(3d0+de)*(a(3,1,i,j+1)+a(3,1,i,j))
     $           *(a(3,2,i,j+1)-a(3,2,i,j))/(df*r)
            fg=fg+(a(1,1,i,j+1)+a(1,1,i,j)
     $           -a(2,2,i,j+1)-a(2,2,i,j))*
     $           (a(2,1,i,j+1)-a(2,1,i,j))/(df*r)
            fg=fg-(a(1,2,i,j+1)+a(1,2,i,j)
     $           +a(2,1,i,j+1)+a(2,1,i,j))*
     $           (a(1,1,i,j+1)-a(1,1,i,j))/(df*r)
            fg=fg+(a(1,3,i,j+1)+a(1,3,i,j))
     $           *(a(2,3,i,j+1)-a(2,3,i,j))/(df*r)
            fg=fg-(a(2,3,i,j+1)+a(2,3,i,j))
     $           *(a(1,3,i,j+1)-a(1,3,i,j))/(df*r)
            fg=fg-(a(3,2,i,j+1)+a(3,2,i,j))
     $           *(a(3,1,i,j+1)-a(3,1,i,j))/(df*r)
 1600    continue
 1700 continue
      do 2000 i=2,n-1
c         do 1900 j=1,m-2
         do 1900 j=2,m-2
            rp=2d0+de
            fg=fg+0.5d0*rp*(a(3,1,i+1,j)-a(3,1,i,j)+a(3,1,i+1,j+1)
     $           -a(3,1,i,j+1))*(a(3,2,i,j+1)-a(3,2,i,j)
     $           +a(3,2,i+1,j+1)-a(3,2,i+1,j))/(dr*df)
            fg=fg+0.5d0*rp*(a(1,1,i+1,j)-a(1,1,i,j)+a(1,1,i+1,j+1)
     $           -a(1,1,i,j+1))*(a(1,2,i,j+1)-a(1,2,i,j)
     $           +a(1,2,i+1,j+1)-a(1,2,i+1,j))/(dr*df)
            fg=fg+0.5d0*rp*(a(2,1,i+1,j)-a(2,1,i,j)+a(2,1,i+1,j+1)
     $           -a(2,1,i,j+1))*(a(2,2,i,j+1)-a(2,2,i,j)
     $           +a(2,2,i+1,j+1)-a(2,2,i+1,j))/(dr*df)
 1900    continue
 2000 continue
      do 2700 i=2,ilimit-1
         j=m-1
         r=dble(i-1)*dr
         do 2500 k=1,3
            do 2400 l=1,3
               if (l.eq.3) then
                  apu=3d0+de
               else
                  apu=1d0
               end if
               fg=fg+(apu/r)*((a(k,l,i,j+1)-a(k,l,i,j))/df)**2
 2400       continue
 2500    continue
         fg=fg+(3d0+de)*(a(1,1,i,j+1)+a(1,1,i,j)
     $        -a(2,2,i,j+1)-a(2,2,i,j))*
     $        (a(1,2,i,j+1)-a(1,2,i,j))/(df*r)
         fg=fg+(3d0+de)*(a(1,2,i,j+1)+a(1,2,i,j)
     $        +a(2,1,i,j+1)+a(2,1,i,j))*
     $        (a(2,2,i,j+1)-a(2,2,i,j))/(df*r)
         fg=fg+(3d0+de)*(a(3,1,i,j+1)+a(3,1,i,j))
     $        *(a(3,2,i,j+1)-a(3,2,i,j))/(df*r)
         fg=fg+(a(1,1,i,j+1)+a(1,1,i,j)
     $        -a(2,2,i,j+1)-a(2,2,i,j))*
     $        (a(2,1,i,j+1)-a(2,1,i,j))/(df*r)
         fg=fg-(a(1,2,i,j+1)+a(1,2,i,j)
     $        +a(2,1,i,j+1)+a(2,1,i,j))*
     $        (a(1,1,i,j+1)-a(1,1,i,j))/(df*r)
         fg=fg+(a(1,3,i,j+1)+a(1,3,i,j))
     $        *(a(2,3,i,j+1)-a(2,3,i,j))/(df*r)
         fg=fg-(a(2,3,i,j+1)+a(2,3,i,j))
     $        *(a(1,3,i,j+1)-a(1,3,i,j))/(df*r)
         fg=fg-(a(3,2,i,j+1)+a(3,2,i,j))
     $        *(a(3,1,i,j+1)-a(3,1,i,j))/(df*r)
 2700 continue
      do 3000 i=2,ilimit-1
         j=m-1
         rp=2d0+de
         fg=fg+0.5d0*rp*(a(3,1,i+1,j)-a(3,1,i,j)+a(3,1,i+1,j+1)
     $        -a(3,1,i,j+1))*(a(3,2,i,j+1)-a(3,2,i,j)
     $        +a(3,2,i+1,j+1)-a(3,2,i+1,j))/(dr*df)
         fg=fg+0.5d0*rp*(a(1,1,i+1,j)-a(1,1,i,j)+a(1,1,i+1,j+1)
     $        -a(1,1,i,j+1))*(a(1,2,i,j+1)-a(1,2,i,j)
     $        +a(1,2,i+1,j+1)-a(1,2,i+1,j))/(dr*df)
         fg=fg+0.5d0*rp*(a(2,1,i+1,j)-a(2,1,i,j)+a(2,1,i+1,j+1)
     $        -a(2,1,i,j+1))*(a(2,2,i,j+1)-a(2,2,i,j)
     $        +a(2,2,i+1,j+1)-a(2,2,i+1,j))/(dr*df)
 3000 continue
      fg=8d0*fg/13d0
      co=24d0/13d0
      do 4000 j=2,m-1
         fsg=fsg-co*0.5d0*(a(1,1,n,j)+a(1,1,n,j))*
     $        (a(1,2,n,j+1)-a(1,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(2,1,n,j)+a(2,1,n,j))*
     $        (a(2,2,n,j+1)-a(2,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(3,1,n,j)+a(3,1,n,j))*
     $        (a(3,2,n,j+1)-a(3,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(1,2,n,j)*a(2,1,n,j)+
     $        a(1,2,n,j+1)*a(2,1,n,j+1))
         fsg=fsg+co*0.5d0*(a(1,1,n,j)*a(2,2,n,j)+
     $        a(1,1,n,j+1)*a(2,2,n,j+1))
 4000 continue
      energy=(fh+fv+fs+fl+fg+fsg)*dr*df
      return
      end

      subroutine gradient(n,m,a,gr)
      implicit none 
      integer i,j,m,n,k,l,nmax,ilimit
      double precision pii,r,dr,df
      parameter (nmax=200,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision gr(3,3,nmax,nmax)
      double precision omega0,de,d0,apu,c,s,co
      double precision ag,rp,rm,l0,ov0,y(nmax),fii
      double precision drr,dfr,drf,dff,drz,dfz,dzr,dzf,dzz
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      ag=8d0/13d0
      ilimit=dint(dble(n)*dsqrt(ov0/omega0)+0.5d0)
      if (omega0.eq.0d0) ilimit=0
      do 40 i=1,n
         do 30 j=1,m
            do 20 k=1,3
               do 10 l=1,3
                  gr(i,j,k,l)=0d0
 10            continue
 20         continue
 30      continue
 40   continue
      rmax=dble(n-1)*dr
      do 200 i=2,n
         r=dr*dble(i-1)
         do 100 j=2,m-1
            if (i.gt.ilimit) then
               v=omega0*r-ov0*rmax*rmax/r
            else
               v=0d0
               gr(3,3,i,j)=gr(3,3,i,j)+2d0*l0*r*a(3,3,i,j)*y(i)
            end if
            gr(3,3,i,j)=gr(3,3,i,j)-4d0*r*y(i)
            gr(3,2,i,j)=gr(3,2,i,j)-4d0*r*y(i)*a(3,2,i,j)*v**2
            gr(1,2,i,j)=gr(1,2,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,2,i,j)+a(2,1,i,j))  
            gr(2,1,i,j)=gr(2,1,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,2,i,j)+a(2,1,i,j))  
            gr(1,1,i,j)=gr(1,1,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,1,i,j)-a(2,2,i,j))  
            gr(2,2,i,j)=gr(2,2,i,j)-(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,1,i,j)-a(2,2,i,j))  
            gr(3,1,i,j)=gr(3,1,i,j)+(2d0*ag*(3d0+de)/r)*y(i)*
     $           a(3,1,i,j)  
            gr(1,3,i,j)=gr(1,3,i,j)+(2d0*ag/r)*a(1,3,i,j)*y(i)  
            gr(3,2,i,j)=gr(3,2,i,j)+(2d0*ag/r)*a(3,2,i,j)*y(i)  
            gr(2,3,i,j)=gr(2,3,i,j)+(2d0*ag/r)*a(2,3,i,j)*y(i)  
 100     continue
 200  continue
      do 300 j=2,m-1
         gr(3,1,n,j)=gr(3,1,n,j)-1d1*(rmax/dr)*d0*a(3,1,n,j)
 300  continue
      do 800 i=2,n-1
         do 700 j=2,m-1
            rp=(dble(i)-0.5d0)*dr
            rm=(dble(i)-1.5d0)*dr
            do 600 k=1,3
               do 500 l=1,3
                  if (l.eq.1) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  gr(k,l,i,j)=gr(k,l,i,j)-2d0*ag*apu*
     $                 (rm*a(k,l,i-1,j)-(rp+rm)*a(k,l,i,j)+
     $                 rp*a(k,l,i+1,j))/(dr**2)
 500           continue
 600        continue
 700     continue
 800  continue
      do 850 i=3,n-1
         do 830 j=2,m-1
            gr(1,1,i,j)=gr(1,1,i,j)+ag*(2d0+de)*
     $           (a(2,2,i+1,j)-a(2,2,i-1,j))/dr
            gr(2,2,i,j)=gr(2,2,i,j)-ag*(2d0+de)*
     $           (a(1,1,i+1,j)-a(1,1,i-1,j))/dr
            gr(1,2,i,j)=gr(1,2,i,j)+ag*(2d0+de)*
     $           (a(2,1,i+1,j)-a(2,1,i-1,j))/dr
            gr(2,1,i,j)=gr(2,1,i,j)-ag*(2d0+de)*
     $           (a(1,2,i+1,j)-a(1,2,i-1,j))/dr
 830     continue
 850  continue
      do 1100 j=2,m-1
         rm=(dble(n)-1.5d0)*dr
         do 1000 k=1,3
            do 900 l=1,3
               if (l.eq.1) then
                  apu=3d0+de
               else
                  apu=1d0
               end if
               gr(k,l,n,j)=gr(k,l,n,j)+2d0*ag*apu*
     $              rm*(a(k,l,n,j)-a(k,l,n-1,j))/(dr**2)
 900        continue
 1000    continue
         gr(1,1,n,j)=gr(1,1,n,j)+ag*(2d0+de)*
     $        (2d0*a(1,1,n,j)-a(2,2,n,j)-a(2,2,n-1,j))/dr
         gr(2,2,n,j)=gr(2,2,n,j)-ag*(2d0+de)*
     $        (a(1,1,n,j)-a(1,1,n-1,j))/dr
         gr(2,1,n,j)=gr(2,1,n,j)+ag*(2d0+de)*
     $        (2d0*a(2,1,n,j)+a(1,2,n,j)+a(1,2,n-1,j))/dr
         gr(1,2,n,j)=gr(1,2,n,j)+ag*(2d0+de)*
     $        (a(2,1,n,j)-a(2,1,n-1,j))/dr
         gr(3,1,n,j)=gr(3,1,n,j)+ag*(2d0+de)*
     $        2d0*a(3,1,n,j)/dr
         gr(1,1,2,j)=gr(1,1,2,j)-ag*(2d0+de)*
     $        (2d0*a(1,1,2,j)-a(2,2,3,j)-a(2,2,2,j))/dr
         gr(2,2,2,j)=gr(2,2,2,j)-ag*(2d0+de)*
     $        (a(1,1,3,j)-a(1,1,2,j))/dr
         gr(2,1,2,j)=gr(2,1,2,j)-ag*(2d0+de)*
     $        (2d0*a(2,1,2,j)+a(1,2,3,j)+a(1,2,2,j))/dr
         gr(1,2,2,j)=gr(1,2,2,j)+ag*(2d0+de)*
     $        (a(2,1,3,j)-a(2,1,2,j))/dr
         gr(3,1,2,j)=gr(3,1,2,j)-ag*(2d0+de)*
     $        2d0*a(3,1,2,j)/dr
 1100 continue
      do 1700 i=2,n
         do 1600 j=2,m-1
            r=dble(i-1)*dr
            do 1500 k=1,3
               do 1400 l=1,3
                  if (l.eq.3) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  gr(k,l,i,j)=gr(k,l,i,j)-2d0*ag*(apu/r)*
     $                 (a(k,l,i,j+1)-2d0*a(k,l,i,j)+
     $                 a(k,l,i,j-1))/(df**2)
 1400          continue
 1500       continue
            gr(1,1,i,j)=gr(1,1,i,j)+(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1))/(df*r)
            gr(2,2,i,j)=gr(2,2,i,j)-(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1))/(df*r)
            gr(1,2,i,j)=gr(1,2,i,j)-(3d0+de)*ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1)-
     $           a(2,2,i,j+1)+a(2,2,i,j-1))/(df*r)
            gr(1,2,i,j)=gr(1,2,i,j)+(3d0+de)*ag*
     $           (a(2,2,i,j+1)-a(2,2,i,j-1))/(df*r)
            gr(2,1,i,j)=gr(2,1,i,j)+(3d0+de)*ag*
     $           (a(2,2,i,j+1)-a(2,2,i,j-1))/(df*r)
            gr(2,2,i,j)=gr(2,2,i,j)-(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+
     $           a(2,1,i,j+1)-a(2,1,i,j-1))/(df*r)
            gr(3,1,i,j)=gr(3,1,i,j)+(3d0+de)*ag*
     $           (a(3,2,i,j+1)-a(3,2,i,j-1))/(df*r)
            gr(3,2,i,j)=gr(3,2,i,j)-(3d0+de)*ag*
     $           (a(3,1,i,j+1)-a(3,1,i,j-1))/(df*r)
            gr(1,1,i,j)=gr(1,1,i,j)+ag*
     $           (a(2,1,i,j+1)-a(2,1,i,j-1))/(df*r)
            gr(2,2,i,j)=gr(2,2,i,j)-ag*
     $           (a(2,1,i,j+1)-a(2,1,i,j-1))/(df*r)
            gr(2,1,i,j)=gr(2,1,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1)-
     $           a(2,2,i,j+1)+a(2,2,i,j-1))/(df*r)
            gr(1,2,i,j)=gr(1,2,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1))/(df*r)
            gr(2,1,i,j)=gr(2,1,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1))/(df*r)
            gr(1,1,i,j)=gr(1,1,i,j)+ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+
     $           a(2,1,i,j+1)-a(2,1,i,j-1))/(df*r)
            gr(1,3,i,j)=gr(1,3,i,j)+ag*
     $           (a(2,3,i,j+1)-a(2,3,i,j-1))/(df*r)
            gr(2,3,i,j)=gr(2,3,i,j)-ag*
     $           (a(1,3,i,j+1)-a(1,3,i,j-1))/(df*r)
            gr(2,3,i,j)=gr(2,3,i,j)-ag*
     $           (a(1,3,i,j+1)-a(1,3,i,j-1))/(df*r)
            gr(1,3,i,j)=gr(1,3,i,j)+ag*
     $           (a(2,3,i,j+1)-a(2,3,i,j-1))/(df*r)
            gr(3,2,i,j)=gr(3,2,i,j)-ag*
     $           (a(3,1,i,j+1)-a(3,1,i,j-1))/(df*r)
            gr(3,1,i,j)=gr(3,1,i,j)+ag*
     $           (a(3,2,i,j+1)-a(3,2,i,j-1))/(df*r)
 1600    continue
 1700 continue
      do 2200 i=3,n-1
         do 2100 j=2,m-1
            gr(3,1,i,j)=gr(3,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(3,2,i,j+1)-a(3,2,i,j-1)+a(3,2,i+1,j+1)
     $           -a(3,2,i+1,j-1))-(a(3,2,i,j+1)-a(3,2,i,j-1)
     $           +a(3,2,i-1,j+1)-a(3,2,i-1,j-1)))/(dr*df)
            gr(3,2,i,j)=gr(3,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(3,1,i+1,j+1)-a(3,1,i+1,j-1)+a(3,1,i,j-1)
     $           -a(3,1,i,j+1))+(a(3,1,i,j+1)-a(3,1,i-1,j+1)
     $           +a(3,1,i-1,j-1)-a(3,1,i,j-1)))/(dr*df)
            gr(1,1,i,j)=gr(1,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+a(1,2,i+1,j+1)
     $           -a(1,2,i+1,j-1))-(a(1,2,i,j+1)-a(1,2,i,j-1)
     $           +a(1,2,i-1,j+1)-a(1,2,i-1,j-1)))/(dr*df)
            gr(1,2,i,j)=gr(1,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(1,1,i+1,j+1)-a(1,1,i+1,j-1)+a(1,1,i,j-1)
     $           -a(1,1,i,j+1))+(a(1,1,i,j+1)-a(1,1,i-1,j+1)
     $           +a(1,1,i-1,j-1)-a(1,1,i,j-1)))/(dr*df)
            gr(2,1,i,j)=gr(2,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(2,2,i,j+1)-a(2,2,i,j-1)+a(2,2,i+1,j+1)
     $           -a(2,2,i+1,j-1))-(a(2,2,i,j+1)-a(2,2,i,j-1)
     $           +a(2,2,i-1,j+1)-a(2,2,i-1,j-1)))/(dr*df)
            gr(2,2,i,j)=gr(2,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(2,1,i+1,j+1)-a(2,1,i+1,j-1)+a(2,1,i,j-1)
     $           -a(2,1,i,j+1))+(a(2,1,i,j+1)-a(2,1,i-1,j+1)
     $           +a(2,1,i-1,j-1)-a(2,1,i,j-1)))/(dr*df)
 2100    continue
 2200 continue
      do 2400 j=2,m-1
         gr(3,1,n,j)=gr(3,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(3,2,n-1,j+1)-a(3,2,n,j+1)+
     $        a(3,2,n-1,j-1)+a(3,2,n,j-1)))/(dr*df)
         gr(3,2,n,j)=gr(3,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(3,1,n,j+1)-a(3,1,n-1,j+1)-
     $        a(3,1,n,j-1)+a(3,1,n-1,j-1)))/(dr*df)
         gr(1,1,n,j)=gr(1,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(1,2,n-1,j+1)-a(1,2,n,j+1)+
     $        a(1,2,n-1,j-1)+a(1,2,n,j-1)))/(dr*df)
         gr(1,2,n,j)=gr(1,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(1,1,n,j+1)-a(1,1,n-1,j+1)-
     $        a(1,1,n,j-1)+a(1,1,n-1,j-1)))/(dr*df)
         gr(2,1,n,j)=gr(2,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(2,2,n-1,j+1)-a(2,2,n,j+1)+
     $        a(2,2,n-1,j-1)+a(2,2,n,j-1)))/(dr*df)
         gr(2,2,n,j)=gr(2,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(2,1,n,j+1)-a(2,1,n-1,j+1)-
     $        a(2,1,n,j-1)+a(2,1,n-1,j-1)))/(dr*df)
         gr(3,1,2,j)=gr(3,1,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(3,2,2,j+1)+a(3,2,3,j+1)-
     $        a(3,2,2,j-1)-a(3,2,3,j-1)))/(dr*df)
         gr(3,2,2,j)=gr(3,2,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(3,1,3,j+1)-a(3,1,3,j-1)-
     $        a(3,1,2,j+1)+a(3,1,2,j-1)))/(dr*df)
         gr(1,1,2,j)=gr(1,1,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(1,2,3,j+1)-a(1,2,3,j-1)+
     $        a(1,2,2,j+1)-a(1,2,2,j-1)))/(dr*df)
         gr(1,2,2,j)=gr(1,2,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(1,1,3,j+1)-a(1,1,3,j-1)-
     $        a(1,1,2,j+1)+a(1,1,2,j-1)))/(dr*df)
         gr(2,1,2,j)=gr(2,1,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(2,2,3,j+1)-a(2,2,3,j-1)+
     $        a(2,2,2,j+1)-a(2,2,2,j-1)))/(dr*df)
         gr(2,2,2,j)=gr(2,2,2,j)-0.5d0*(2d0+de)*ag*(
     $        (a(2,1,3,j+1)-a(2,1,3,j-1)-
     $        a(2,1,2,j+1)+a(2,1,2,j-1)))/(dr*df)
 2400 continue
      do 3000 j=2,m-1
         apu=-2d0*ag*0.5d0*df
         fii=(dble(j)-1.5d0)*df
         c=cos(fii)
         s=sin(fii)
         drr=apu*(3d0+de)*(a(1,1,2,j)-a(1,1,1,j))
         drf=apu*(a(1,2,2,j)-a(1,2,1,j))
         dfr=apu*(3d0+de)*(a(2,1,2,j)-a(2,1,1,j))
         dff=apu*(a(2,2,2,j)-a(2,2,1,j))
         drz=apu*(a(1,3,2,j)-a(1,3,1,j))
         dzr=apu*(3d0+de)*(a(3,1,2,j)-a(3,1,1,j))
         dfz=apu*(a(2,3,2,j)-a(2,3,1,j))
         dzf=apu*(a(3,2,2,j)-a(3,2,1,j))
         dzz=apu*(a(3,3,2,j)-a(3,3,1,j))
         gr(1,1,1,1)=gr(1,1,1,1)+(c*c*drr-s*c*(dfr+drf)+s*s*dff)
         gr(1,2,1,1)=gr(1,2,1,1)+(c*c*drf+s*c*(drr-dff)-s*s*dfr)
         gr(1,3,1,1)=gr(1,3,1,1)+(c*drz-s*dfz)
         gr(2,1,1,1)=gr(2,1,1,1)+(c*c*dfr+s*c*(drr-dff)-s*s*drf)
         gr(2,2,1,1)=gr(2,2,1,1)+(c*c*dff+s*c*(dfr+drf)+s*s*drr)
         gr(2,3,1,1)=gr(2,3,1,1)+(c*dfz+s*drz)
         gr(3,1,1,1)=gr(3,1,1,1)+(c*dzr-s*dzf)
         gr(3,2,1,1)=gr(3,2,1,1)+(c*dzf+s*dzr)
         gr(3,3,1,1)=gr(3,3,1,1)+dzz
 3000 continue
      co=24d0/13d0
      do 3001 j=2,m-1
         gr(1,1,n,j)=gr(1,1,n,j)-0.5d0*co*
     $        (a(1,2,n,j+1)-a(1,2,n,j-1))/df
         gr(1,2,n,j)=gr(1,2,n,j)-0.5d0*co*
     $        (a(1,1,n,j+1)-a(1,1,n,j-1))/df
         gr(2,1,n,j)=gr(2,1,n,j)-0.5d0*co*
     $        (a(2,2,n,j+1)-a(2,2,n,j-1))/df
         gr(2,2,n,j)=gr(2,2,n,j)-0.5d0*co*
     $        (a(2,1,n,j+1)-a(2,1,n,j-1))/df
         gr(3,1,n,j)=gr(3,1,n,j)-0.5d0*co*
     $        (a(3,2,n,j+1)-a(3,2,n,j-1))/df
         gr(3,2,n,j)=gr(3,2,n,j)-0.5d0*co*
     $        (a(3,1,n,j+1)-a(3,1,n,j-1))/df
         gr(1,2,n,j)=gr(1,2,n,j)-co*a(2,1,n,j)
         gr(2,1,n,j)=gr(2,1,n,j)-co*a(1,2,n,j)
         gr(1,1,n,j)=gr(1,1,n,j)+co*a(2,2,n,j)
         gr(2,2,n,j)=gr(2,2,n,j)+co*a(1,1,n,j)
 3001 continue
c      write (6,*) drr
c      write (6,*) drf
c      write (6,*) dfr
c      write (6,*) dff
c      write (6,*) drz
c      write (6,*) dzr
c      write (6,*) dfz
c      write (6,*) dzf
c      write (6,*) dzz
      do 3040 i=2,n
         do 3030 j=2,m-1
            do 3020 k=1,3
               do 3010 l=1,3
                  gr(k,l,i,j)=gr(k,l,i,j)*dr*df
 3010          continue
 3020       continue
 3030    continue
 3040 continue
      do 3070 k=1,3
         do 3060 l=1,3
            gr(k,l,1,1)=gr(k,l,1,1)*dr*df
 3060    continue
 3070 continue
      do 4040 i=ilimit,n
         do 4020 k=1,3
            do 4010 l=1,3
               gr(k,l,i,2)=0d0
               gr(k,l,i,m-1)=0d0
 4010       continue
 4020    continue
 4040 continue
      return
      end

      subroutine energy2(n,m,a)
      implicit none 
      integer i,j,m,n,k,l,nmax,ilimit
      double precision pii,r,dr,df,rp
      parameter (nmax=200,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision omega0,de,d0,apu,eps,fl,fsg,co
      double precision fh,fv,fs,fg,l0,ov0,y(nmax)
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      fh=0d0
      fv=0d0
      fg=0d0
      fs=0d0
      fl=0d0
      fsg=0d0
      eps=1d-8
      rmax=dble(n-1)*dr
      ilimit=dint(dble(n)*dsqrt(ov0/omega0)+0.5d0)
      if (omega0.eq.0d0) ilimit=0
      do 200 i=2,n
         r=dr*dble(i-1)
         do 100 j=2,m-1
            if (i.gt.ilimit) then
               v=omega0*r-ov0*rmax*rmax/r
            else
               v=0d0
               fl=fl+l0*r*y(i)*a(3,3,i,j)**2
            end if
            fh=fh-r*(4d0*a(3,3,i,j)+1d0)*y(i)
            fv=fv-2d0*r*y(i)*(a(3,2,i,j)*v)**2
            fg=fg+(1d0/r)*y(i)*(4d0+de)*
     $           ((a(1,2,i,j)+a(2,1,i,j))**2+
     $           (a(1,1,i,j)-a(2,2,i,j))**2)
            fg=fg+(1d0/r)*y(i)*((3d0+de)*
     $           a(3,1,i,j)**2+a(1,3,i,j)**2+
     $           a(3,2,i,j)**2+a(2,3,i,j)**2)
 100     continue
 200  continue
      do 300 j=2,m-1
         fs=fs-5d0*(rmax/dr)*d0*a(3,1,n,j)**2
 300  continue
      do 800 i=1,n-1
         do 700 j=2,m-1
            rp=(dble(i)-0.5d0)*dr
            do 600 k=1,3
               do 500 l=1,3
                  if (l.eq.1) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  fg=fg+apu*rp*((a(k,l,i+1,j)-a(k,l,i,j))/dr)**2
 500           continue
 600        continue
            fg=fg+(2d0+de)*(a(1,1,i+1,j)+a(1,1,i,j)
     $           -a(2,2,i+1,j)-a(2,2,i,j))*
     $           (a(1,1,i+1,j)-a(1,1,i,j))/dr
            fg=fg+(2d0+de)*(a(1,2,i+1,j)+a(1,2,i,j)
     $           +a(2,1,i+1,j)+a(2,1,i,j))*
     $           (a(2,1,i+1,j)-a(2,1,i,j))/dr
            fg=fg+(2d0+de)*(a(3,1,i+1,j)+a(3,1,i,j))
     $           *(a(3,1,i+1,j)-a(3,1,i,j))/dr
 700     continue
 800  continue
      do 1700 i=2,n
         do 1600 j=1,m-2
c         do 1600 j=2,m-2
            r=dble(i-1)*dr
            do 1500 k=1,3
               do 1400 l=1,3
                  if (l.eq.3) then
                     apu=3d0+de
                  else
                     apu=1d0
                  end if
                  fg=fg+(apu/r)*((a(k,l,i,j+1)-a(k,l,i,j))/df)**2
 1400          continue
 1500       continue
            fg=fg+(3d0+de)*(a(1,1,i,j+1)+a(1,1,i,j)
     $           -a(2,2,i,j+1)-a(2,2,i,j))*
     $           (a(1,2,i,j+1)-a(1,2,i,j))/(df*r)
            fg=fg+(3d0+de)*(a(1,2,i,j+1)+a(1,2,i,j)
     $           +a(2,1,i,j+1)+a(2,1,i,j))*
     $           (a(2,2,i,j+1)-a(2,2,i,j))/(df*r)
            fg=fg+(3d0+de)*(a(3,1,i,j+1)+a(3,1,i,j))
     $           *(a(3,2,i,j+1)-a(3,2,i,j))/(df*r)
            fg=fg+(a(1,1,i,j+1)+a(1,1,i,j)
     $           -a(2,2,i,j+1)-a(2,2,i,j))*
     $           (a(2,1,i,j+1)-a(2,1,i,j))/(df*r)
            fg=fg-(a(1,2,i,j+1)+a(1,2,i,j)
     $           +a(2,1,i,j+1)+a(2,1,i,j))*
     $           (a(1,1,i,j+1)-a(1,1,i,j))/(df*r)
            fg=fg+(a(1,3,i,j+1)+a(1,3,i,j))
     $           *(a(2,3,i,j+1)-a(2,3,i,j))/(df*r)
            fg=fg-(a(2,3,i,j+1)+a(2,3,i,j))
     $           *(a(1,3,i,j+1)-a(1,3,i,j))/(df*r)
            fg=fg-(a(3,2,i,j+1)+a(3,2,i,j))
     $           *(a(3,1,i,j+1)-a(3,1,i,j))/(df*r)
 1600    continue
 1700 continue
      do 2000 i=1,n-1
         do 1900 j=1,m-2
c         do 1900 j=2,m-2
            rp=2d0+de
            fg=fg+0.5d0*rp*(a(3,1,i+1,j)-a(3,1,i,j)+a(3,1,i+1,j+1)
     $           -a(3,1,i,j+1))*(a(3,2,i,j+1)-a(3,2,i,j)
     $           +a(3,2,i+1,j+1)-a(3,2,i+1,j))/(dr*df)
            fg=fg+0.5d0*rp*(a(1,1,i+1,j)-a(1,1,i,j)+a(1,1,i+1,j+1)
     $           -a(1,1,i,j+1))*(a(1,2,i,j+1)-a(1,2,i,j)
     $           +a(1,2,i+1,j+1)-a(1,2,i+1,j))/(dr*df)
            fg=fg+0.5d0*rp*(a(2,1,i+1,j)-a(2,1,i,j)+a(2,1,i+1,j+1)
     $           -a(2,1,i,j+1))*(a(2,2,i,j+1)-a(2,2,i,j)
     $           +a(2,2,i+1,j+1)-a(2,2,i+1,j))/(dr*df)
 1900    continue
 2000 continue
      do 2700 i=2,ilimit-1
         j=m-1
         r=dble(i-1)*dr
         do 2500 k=1,3
            do 2400 l=1,3
               if (l.eq.3) then
                  apu=3d0+de
               else
                  apu=1d0
               end if
               fg=fg+(apu/r)*((a(k,l,i,j+1)-a(k,l,i,j))/df)**2
 2400       continue
 2500    continue
         fg=fg+(3d0+de)*(a(1,1,i,j+1)+a(1,1,i,j)
     $        -a(2,2,i,j+1)-a(2,2,i,j))*
     $        (a(1,2,i,j+1)-a(1,2,i,j))/(df*r)
         fg=fg+(3d0+de)*(a(1,2,i,j+1)+a(1,2,i,j)
     $        +a(2,1,i,j+1)+a(2,1,i,j))*
     $        (a(2,2,i,j+1)-a(2,2,i,j))/(df*r)
         fg=fg+(3d0+de)*(a(3,1,i,j+1)+a(3,1,i,j))
     $        *(a(3,2,i,j+1)-a(3,2,i,j))/(df*r)
         fg=fg+(a(1,1,i,j+1)+a(1,1,i,j)
     $        -a(2,2,i,j+1)-a(2,2,i,j))*
     $        (a(2,1,i,j+1)-a(2,1,i,j))/(df*r)
         fg=fg-(a(1,2,i,j+1)+a(1,2,i,j)
     $        +a(2,1,i,j+1)+a(2,1,i,j))*
     $        (a(1,1,i,j+1)-a(1,1,i,j))/(df*r)
         fg=fg+(a(1,3,i,j+1)+a(1,3,i,j))
     $        *(a(2,3,i,j+1)-a(2,3,i,j))/(df*r)
         fg=fg-(a(2,3,i,j+1)+a(2,3,i,j))
     $        *(a(1,3,i,j+1)-a(1,3,i,j))/(df*r)
         fg=fg-(a(3,2,i,j+1)+a(3,2,i,j))
     $        *(a(3,1,i,j+1)-a(3,1,i,j))/(df*r)
 2700 continue
      do 3000 i=1,ilimit-1
         j=m-1
         rp=2d0+de
         fg=fg+0.5d0*rp*(a(3,1,i+1,j)-a(3,1,i,j)+a(3,1,i+1,j+1)
     $        -a(3,1,i,j+1))*(a(3,2,i,j+1)-a(3,2,i,j)
     $        +a(3,2,i+1,j+1)-a(3,2,i+1,j))/(dr*df)
         fg=fg+0.5d0*rp*(a(1,1,i+1,j)-a(1,1,i,j)+a(1,1,i+1,j+1)
     $        -a(1,1,i,j+1))*(a(1,2,i,j+1)-a(1,2,i,j)
     $        +a(1,2,i+1,j+1)-a(1,2,i+1,j))/(dr*df)
         fg=fg+0.5d0*rp*(a(2,1,i+1,j)-a(2,1,i,j)+a(2,1,i+1,j+1)
     $        -a(2,1,i,j+1))*(a(2,2,i,j+1)-a(2,2,i,j)
     $        +a(2,2,i+1,j+1)-a(2,2,i+1,j))/(dr*df)
 3000 continue
      fg=8d0*fg/13d0
      co=24d0/13d0
      do 4000 j=2,m-1
         fsg=fsg-co*0.5d0*(a(1,1,n,j)+a(1,1,n,j))*
     $        (a(1,2,n,j+1)-a(1,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(2,1,n,j)+a(2,1,n,j))*
     $        (a(2,2,n,j+1)-a(2,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(3,1,n,j)+a(3,1,n,j))*
     $        (a(3,2,n,j+1)-a(3,2,n,j))/df
         fsg=fsg-co*0.5d0*(a(1,2,n,j)*a(2,1,n,j)+
     $        a(1,2,n,j+1)*a(2,1,n,j+1))
         fsg=fsg+co*0.5d0*(a(1,1,n,j)*a(2,2,n,j)+
     $        a(1,1,n,j+1)*a(2,2,n,j+1))
 4000 continue      
      write (6,*) 'Magnetic energy =',fh*dr*df
      write (6,*) 'Flow energy =',fv*dr*df
      write (6,*) 'Surface energy =',fs*dr*df
      write (6,*) 'Gradient energy =',fg*dr*df
      write (6,*) 'Vortex energy =',fl*dr*df
      write (6,*) 'New energy =',fsg*df
      return
      end
