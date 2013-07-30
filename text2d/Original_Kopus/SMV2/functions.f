
      double precision function energy(n,m,a)
      implicit none 
      integer i,j,m,n,k,l,nmax,ilimit
      double precision pii,r,dr,df,rp
      parameter (nmax=500,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision omega0,de,d0,apu,eps
      double precision fh,fv,fs,fg,l0,fl,ov0,y(nmax)
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      fh=0d0
      fv=0d0
      fg=0d0
      fs=0d0
      fl=0d0
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
         do 1600 j=1,m-1
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
         do 1900 j=1,m-1
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
      fg=8d0*fg/13d0
      energy=(fh+fv+fs+fl+fg)*dr*df
      return
      end

      subroutine gradient(n,m,a,gr)
      implicit none 
      integer i,j,m,n,k,l,nmax,ilimit
      double precision pii,r,dr,df
      parameter (nmax=500,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision gr(3,3,nmax,nmax)
      double precision omega0,de,d0,apu
      double precision ag,rp,rm,l0,ov0,y(nmax)
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
               gr(3,3,i,j)=gr(3,3,i,j)+2d0*l0*r*a(3,3,i,j)*y(i)*dr*df
            end if
            gr(3,3,i,j)=gr(3,3,i,j)-4d0*r*y(i)*dr*df
            gr(3,2,i,j)=gr(3,2,i,j)-4d0*r*y(i)*a(3,2,i,j)*v**2*dr*df
            gr(1,2,i,j)=gr(1,2,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,2,i,j)+a(2,1,i,j))*dr*df  
            gr(2,1,i,j)=gr(2,1,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,2,i,j)+a(2,1,i,j))*dr*df  
            gr(1,1,i,j)=gr(1,1,i,j)+(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,1,i,j)-a(2,2,i,j))*dr*df  
            gr(2,2,i,j)=gr(2,2,i,j)-(2d0*ag*(4d0+de)/r)*y(i)*
     $           (a(1,1,i,j)-a(2,2,i,j))*dr*df  
            gr(3,1,i,j)=gr(3,1,i,j)+(2d0*ag*(3d0+de)/r)*y(i)*
     $           a(3,1,i,j)*dr*df 
            gr(1,3,i,j)=gr(1,3,i,j)+(2d0*ag/r)*a(1,3,i,j)*y(i)*dr*df  
            gr(3,2,i,j)=gr(3,2,i,j)+(2d0*ag/r)*a(3,2,i,j)*y(i)*dr*df  
            gr(2,3,i,j)=gr(2,3,i,j)+(2d0*ag/r)*a(2,3,i,j)*y(i)*dr*df  
 100     continue
 200  continue
      do 300 j=2,m-1
         gr(3,1,n,j)=gr(3,1,n,j)-1d1*rmax*d0*a(3,1,n,j)*df
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
     $                 rp*a(k,l,i+1,j))*df/dr
 500           continue
 600        continue
            gr(1,1,i,j)=gr(1,1,i,j)+ag*(2d0+de)*
     $           (a(2,2,i+1,j)-a(2,2,i-1,j))*df
            gr(2,2,i,j)=gr(2,2,i,j)-ag*(2d0+de)*
     $           (a(1,1,i+1,j)-a(1,1,i-1,j))*df
            gr(1,2,i,j)=gr(1,2,i,j)+ag*(2d0+de)*
     $           (a(2,1,i+1,j)-a(2,1,i-1,j))*df
            gr(2,1,i,j)=gr(2,1,i,j)-ag*(2d0+de)*
     $           (a(1,2,i+1,j)-a(1,2,i-1,j))*df
 700     continue
 800  continue
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
     $              rm*(a(k,l,n,j)-a(k,l,n-1,j))*df/dr
 900        continue
 1000    continue
         gr(1,1,n,j)=gr(1,1,n,j)+ag*(2d0+de)*
     $        (2d0*a(1,1,n,j)-a(2,2,n,j)-a(2,2,n-1,j))*df
         gr(2,2,n,j)=gr(2,2,n,j)-ag*(2d0+de)*
     $        (a(1,1,n,j)-a(1,1,n-1,j))*df
         gr(2,1,n,j)=gr(2,1,n,j)+ag*(2d0+de)*
     $        (2d0*a(2,1,n,j)+a(1,2,n,j)+a(1,2,n-1,j))*df
         gr(1,2,n,j)=gr(1,2,n,j)+ag*(2d0+de)*
     $        (a(2,1,n,j)-a(2,1,n-1,j))*df
         gr(3,1,n,j)=gr(3,1,n,j)+ag*(2d0+de)*
     $        2d0*a(3,1,n,j)*df
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
     $                 a(k,l,i,j-1))*dr/df
 1400          continue
 1500       continue
            gr(1,1,i,j)=gr(1,1,i,j)+(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1))*dr/r
            gr(2,2,i,j)=gr(2,2,i,j)-(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1))*dr/r
            gr(1,2,i,j)=gr(1,2,i,j)-(3d0+de)*ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1)-
     $           a(2,2,i,j+1)+a(2,2,i,j-1))*dr/r
            gr(1,2,i,j)=gr(1,2,i,j)+(3d0+de)*ag*
     $           (a(2,2,i,j+1)-a(2,2,i,j-1))*dr/r
            gr(2,1,i,j)=gr(2,1,i,j)+(3d0+de)*ag*
     $           (a(2,2,i,j+1)-a(2,2,i,j-1))*dr/r
            gr(2,2,i,j)=gr(2,2,i,j)-(3d0+de)*ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+
     $           a(2,1,i,j+1)-a(2,1,i,j-1))*dr/r
            gr(3,1,i,j)=gr(3,1,i,j)+(3d0+de)*ag*
     $           (a(3,2,i,j+1)-a(3,2,i,j-1))*dr/r
            gr(3,2,i,j)=gr(3,2,i,j)-(3d0+de)*ag*
     $           (a(3,1,i,j+1)-a(3,1,i,j-1))*dr/r
            gr(1,1,i,j)=gr(1,1,i,j)+ag*
     $           (a(2,1,i,j+1)-a(2,1,i,j-1))*dr/r
            gr(2,2,i,j)=gr(2,2,i,j)-ag*
     $           (a(2,1,i,j+1)-a(2,1,i,j-1))*dr/r
            gr(2,1,i,j)=gr(2,1,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1)-
     $           a(2,2,i,j+1)+a(2,2,i,j-1))*dr/r
            gr(1,2,i,j)=gr(1,2,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1))*dr/r
            gr(2,1,i,j)=gr(2,1,i,j)-ag*
     $           (a(1,1,i,j+1)-a(1,1,i,j-1))*dr/r
            gr(1,1,i,j)=gr(1,1,i,j)+ag*
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+
     $           a(2,1,i,j+1)-a(2,1,i,j-1))*dr/r
            gr(1,3,i,j)=gr(1,3,i,j)+ag*
     $           (a(2,3,i,j+1)-a(2,3,i,j-1))*dr/r
            gr(2,3,i,j)=gr(2,3,i,j)-ag*
     $           (a(1,3,i,j+1)-a(1,3,i,j-1))*dr/r
            gr(2,3,i,j)=gr(2,3,i,j)-ag*
     $           (a(1,3,i,j+1)-a(1,3,i,j-1))*dr/r
            gr(1,3,i,j)=gr(1,3,i,j)+ag*
     $           (a(2,3,i,j+1)-a(2,3,i,j-1))*dr/r
            gr(3,2,i,j)=gr(3,2,i,j)-ag*
     $           (a(3,1,i,j+1)-a(3,1,i,j-1))*dr/r
            gr(3,1,i,j)=gr(3,1,i,j)+ag*
     $           (a(3,2,i,j+1)-a(3,2,i,j-1))*dr/r
 1600    continue
 1700 continue
      do 2200 i=2,n-1
         do 2100 j=2,m-1
            gr(3,1,i,j)=gr(3,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(3,2,i,j+1)-a(3,2,i,j-1)+a(3,2,i+1,j+1)
     $           -a(3,2,i+1,j-1))-(a(3,2,i,j+1)-a(3,2,i,j-1)
     $           +a(3,2,i-1,j+1)-a(3,2,i-1,j-1)))
            gr(3,2,i,j)=gr(3,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(3,1,i+1,j+1)-a(3,1,i+1,j-1)+a(3,1,i,j-1)
     $           -a(3,1,i,j+1))+(a(3,1,i,j+1)-a(3,1,i-1,j+1)
     $           +a(3,1,i-1,j-1)-a(3,1,i,j-1)))
            gr(1,1,i,j)=gr(1,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(1,2,i,j+1)-a(1,2,i,j-1)+a(1,2,i+1,j+1)
     $           -a(1,2,i+1,j-1))-(a(1,2,i,j+1)-a(1,2,i,j-1)
     $           +a(1,2,i-1,j+1)-a(1,2,i-1,j-1)))
            gr(1,2,i,j)=gr(1,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(1,1,i+1,j+1)-a(1,1,i+1,j-1)+a(1,1,i,j-1)
     $           -a(1,1,i,j+1))+(a(1,1,i,j+1)-a(1,1,i-1,j+1)
     $           +a(1,1,i-1,j-1)-a(1,1,i,j-1)))
            gr(2,1,i,j)=gr(2,1,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(2,2,i,j+1)-a(2,2,i,j-1)+a(2,2,i+1,j+1)
     $           -a(2,2,i+1,j-1))-(a(2,2,i,j+1)-a(2,2,i,j-1)
     $           +a(2,2,i-1,j+1)-a(2,2,i-1,j-1)))
            gr(2,2,i,j)=gr(2,2,i,j)-0.5d0*(2d0+de)*ag*(
     $           (a(2,1,i+1,j+1)-a(2,1,i+1,j-1)+a(2,1,i,j-1)
     $           -a(2,1,i,j+1))+(a(2,1,i,j+1)-a(2,1,i-1,j+1)
     $           +a(2,1,i-1,j-1)-a(2,1,i,j-1)))
 2100    continue
 2200 continue
c$$$      do 2300 i=2,n-1
c$$$         gr(3,1,i,1)=gr(3,1,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(3,2,i,2)+a(3,2,i,1)-a(3,2,i+1,2)+a(3,2,i+1,1))
c$$$     $        +(a(3,2,i-1,2)-a(3,2,i-1,1)+a(3,2,i,2)
c$$$     $        -a(3,2,i,1)))/(dr*df)
c$$$         gr(3,1,i,m)=gr(3,1,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(3,2,i,m)+a(3,2,i,m-1)-a(3,2,i+1,m)+a(3,2,i+1,m-1))
c$$$     $        +(a(3,2,i-1,m)-a(3,2,i-1,m-1)+a(3,2,i,m)
c$$$     $        -a(3,2,i,m-1)))/(dr*df)
c$$$         gr(3,2,i,1)=gr(3,2,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(3,1,i+1,1)+a(3,1,i,1)-a(3,1,i+1,2)+a(3,1,i,2))
c$$$     $        -(a(3,1,i,1)-a(3,1,i-1,1)+a(3,1,i,2)
c$$$     $        -a(3,1,i-1,2)))/(dr*df)
c$$$         gr(3,2,i,m)=gr(3,2,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (a(3,1,i+1,m-1)-a(3,1,i,m-1)+a(3,1,i+1,m)-a(3,1,i,m))
c$$$     $        +(a(3,1,i,m-1)-a(3,1,i-1,m-1)+a(3,1,i,m)
c$$$     $        -a(3,1,i-1,m)))/(dr*df)
c$$$         gr(1,1,i,1)=gr(1,1,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(1,2,i,2)+a(1,2,i,1)-a(1,2,i+1,2)+a(1,2,i+1,1))
c$$$     $        +(a(1,2,i-1,2)-a(1,2,i-1,1)+a(1,2,i,2)
c$$$     $        -a(1,2,i,1)))/(dr*df)
c$$$         gr(1,1,i,m)=gr(1,1,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(1,2,i,m)+a(1,2,i,m-1)-a(1,2,i+1,m)+a(1,2,i+1,m-1))
c$$$     $        +(a(1,2,i-1,m)-a(1,2,i-1,m-1)+a(1,2,i,m)
c$$$     $        -a(1,2,i,m-1)))/(dr*df)
c$$$         gr(1,2,i,1)=gr(1,2,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(1,1,i+1,1)+a(1,1,i,1)-a(1,1,i+1,2)+a(1,1,i,2))
c$$$     $        -(a(1,1,i,1)-a(1,1,i-1,1)+a(1,1,i,2)
c$$$     $        -a(1,1,i-1,2)))/(dr*df)
c$$$         gr(1,2,i,m)=gr(1,2,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (a(1,1,i+1,m-1)-a(1,1,i,m-1)+a(1,1,i+1,m)-a(1,1,i,m))
c$$$     $        +(a(1,1,i,m-1)-a(1,1,i-1,m-1)+a(1,1,i,m)
c$$$     $        -a(1,1,i-1,m)))/(dr*df)
c$$$         gr(2,1,i,1)=gr(2,1,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(2,2,i,2)+a(2,2,i,1)-a(2,2,i+1,2)+a(2,2,i+1,1))
c$$$     $        +(a(2,2,i-1,2)-a(2,2,i-1,1)+a(2,2,i,2)
c$$$     $        -a(2,2,i,1)))/(dr*df)
c$$$         gr(2,1,i,m)=gr(2,1,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(2,2,i,m)+a(2,2,i,m-1)-a(2,2,i+1,m)+a(2,2,i+1,m-1))
c$$$     $        +(a(2,2,i-1,m)-a(2,2,i-1,m-1)+a(2,2,i,m)
c$$$     $        -a(2,2,i,m-1)))/(dr*df)
c$$$         gr(2,2,i,1)=gr(2,2,i,1)+0.5d0*(2d0+de)*ag*(
c$$$     $        (-a(2,1,i+1,1)+a(2,1,i,1)-a(2,1,i+1,2)+a(2,1,i,2))
c$$$     $        -(a(2,1,i,1)-a(2,1,i-1,1)+a(2,1,i,2)
c$$$     $        -a(2,1,i-1,2)))/(dr*df)
c$$$         gr(2,2,i,m)=gr(2,2,i,m)+0.5d0*(2d0+de)*ag*(
c$$$     $        (a(2,1,i+1,m-1)-a(2,1,i,m-1)+a(2,1,i+1,m)-a(2,1,i,m))
c$$$     $        +(a(2,1,i,m-1)-a(2,1,i-1,m-1)+a(2,1,i,m)
c$$$     $        -a(2,1,i-1,m)))/(dr*df)
c$$$ 2300 continue
      do 2400 j=2,m-1
         gr(3,1,n,j)=gr(3,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(3,2,n-1,j+1)-a(3,2,n,j+1)+
     $        a(3,2,n-1,j-1)+a(3,2,n,j-1)))
         gr(3,2,n,j)=gr(3,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(3,1,n,j+1)-a(3,1,n-1,j+1)-
     $        a(3,1,n,j-1)+a(3,1,n-1,j-1)))
         gr(1,1,n,j)=gr(1,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(1,2,n-1,j+1)-a(1,2,n,j+1)+
     $        a(1,2,n-1,j-1)+a(1,2,n,j-1)))
         gr(1,2,n,j)=gr(1,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(1,1,n,j+1)-a(1,1,n-1,j+1)-
     $        a(1,1,n,j-1)+a(1,1,n-1,j-1)))
         gr(2,1,n,j)=gr(2,1,n,j)-0.5d0*(2d0+de)*ag*(
     $        (-a(2,2,n-1,j+1)-a(2,2,n,j+1)+
     $        a(2,2,n-1,j-1)+a(2,2,n,j-1)))
         gr(2,2,n,j)=gr(2,2,n,j)-0.5d0*(2d0+de)*ag*(
     $        (a(2,1,n,j+1)-a(2,1,n-1,j+1)-
     $        a(2,1,n,j-1)+a(2,1,n-1,j-1)))
 2400 continue
c$$$      gr(3,1,n,1)=gr(3,1,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (a(3,2,n-1,2)-a(3,2,n-1,1)+
c$$$     $     a(3,2,n,2)-a(3,2,n,1)))/(dr*df)
c$$$      gr(3,2,n,1)=gr(3,2,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(3,1,n,1)+a(3,1,n-1,1)-
c$$$     $     a(3,1,n,2)+a(3,1,n-1,2)))/(dr*df)
c$$$      gr(3,1,n,m)=gr(3,1,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(3,2,n-1,m)+a(3,2,n-1,m-1)+
c$$$     $     -a(3,2,n,m)+a(3,2,n,m-1)))/(dr*df)
c$$$      gr(3,2,n,m)=gr(3,2,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(3,1,n,m-1)+a(3,1,n-1,m-1)-
c$$$     $     a(3,1,n,m)+a(3,1,n-1,m)))/(dr*df)
c$$$      gr(1,1,n,1)=gr(1,1,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (a(1,2,n-1,2)-a(1,2,n-1,1)+
c$$$     $     a(1,2,n,2)-a(1,2,n,1)))/(dr*df)
c$$$      gr(1,2,n,1)=gr(1,2,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(1,1,n,1)+a(1,1,n-1,1)-
c$$$     $     a(1,1,n,2)+a(1,1,n-1,2)))/(dr*df)
c$$$      gr(1,1,n,m)=gr(1,1,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(1,2,n-1,m)+a(1,2,n-1,m-1)+
c$$$     $     -a(1,2,n,m)+a(1,2,n,m-1)))/(dr*df)
c$$$      gr(1,2,n,m)=gr(1,2,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(1,1,n,m-1)+a(1,1,n-1,m-1)-
c$$$     $     a(1,1,n,m)+a(1,1,n-1,m)))/(dr*df)
c$$$      gr(2,1,n,1)=gr(2,1,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (a(2,2,n-1,2)-a(2,2,n-1,1)+
c$$$     $     a(2,2,n,2)-a(2,2,n,1)))/(dr*df)
c$$$      gr(2,2,n,1)=gr(2,2,n,1)+0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(2,1,n,1)+a(2,1,n-1,1)-
c$$$     $     a(2,1,n,2)+a(2,1,n-1,2)))/(dr*df)
c$$$      gr(2,1,n,m)=gr(2,1,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(2,2,n-1,m)+a(2,2,n-1,m-1)+
c$$$     $     -a(2,2,n,m)+a(2,2,n,m-1)))/(dr*df)
c$$$      gr(2,2,n,m)=gr(2,2,n,m)-0.5d0*(2d0+de)*ag*(
c$$$     $     (-a(2,1,n,m-1)+a(2,1,n-1,m-1)-
c$$$     $     a(2,1,n,m)+a(2,1,n-1,m)))/(dr*df)
c      do 3040 i=2,n
c         do 3030 j=2,m-1
c            do 3020 k=1,3
c               do 3010 l=1,3
c                  gr(k,l,i,j)=gr(k,l,i,j)*dr*df
c 3010          continue
c 3020       continue
c 3030    continue
c 3040 continue
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
      parameter (nmax=500,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),v,rmax
      double precision omega0,de,d0,apu,eps,fl
      double precision fh,fv,fs,fg,l0,ov0,y(nmax)
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      fh=0d0
      fv=0d0
      fg=0d0
      fs=0d0
      fl=0d0
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
         do 1600 j=1,m-1
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
         do 1900 j=1,m-1
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
      fg=8d0*fg/13d0
      write (6,*) 'Magnetic energy =',fh*dr*df
      write (6,*) 'Flow energy =',fv*dr*df
      write (6,*) 'Surface energy =',fs*dr*df
      write (6,*) 'Gradient energy =',fg*dr*df
      write (6,*) 'Vortex energy =',fl*dr*df
      return
      end
