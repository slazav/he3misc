      program Bphase2D
      implicit none 
      integer i,j,m,n,nmax,query,iter,it,number,k,l,iterm,itm,ilimit,jj
      integer ii,kk,napu
      double precision pii,dr,df,vaim1,vaim2,summa
      parameter (n=20,m=75,number=100)
      parameter (nmax=200,pii=3.14159265359d0)
      double precision a(3,3,nmax,nmax),gr(3,3,nmax,nmax)
      double precision nvf(3,nmax,nmax),nlen,dnx,dny,dnz,ksumma(m)
      double precision t(3,nmax,nmax),big,nmr(number+1)
      double precision rmax,len,beta,alpha,nv(3,nmax,nmax)
      double precision omega,de,dpera,eps,energy,lp,nyyb,dhh,tau
      double precision vd,ksiih,e2,e1,temp,p,h,sum,iso,fii,apu
      double precision con(-2*number:2*number),ov,vor,d0,y(nmax)
      double precision omega0,ov0,l0,xx,yy,japu
      double precision kcon(m,-2*number:2*number),knmr(m,number+1)
      common /param/ dr,df,d0,omega0,de,l0,ov0,y      
      open (11,file='alkuarvot.dat')
      read (11,*) query
      read (11,*) eps
      read (11,*) rmax
      read (11,*) p
      read (11,*) omega
      read (11,*) vor
      read (11,*) h
      read (11,*) temp
      read (11,*) nyyb
      read (11,*) dhh
      read (11,*) tau
      read (11,*) iterm
      read (11,*) itm
      do 10 i=2,n
         y(i)=1d0
         if (i.eq.n) y(i)=0.5d0
 10   continue
      ov=1.05d-4*dble(vor)/(rmax*rmax)
      call vakiot(temp,omega,p,h,rmax,dpera,vd,lp,ksiih,de,dhh,tau,
     $     nyyb,vaim1,vaim2)
      write (6,*) 'Vd (mm/s)=',vd*10d0
      write (6,*) 'Delta =',de
      write (6,*) 'd/a (mm)=',dpera*10d0
      write (6,*) 'Xi h (mm)=',ksiih*10d0
      write (6,*) 'lp =',lp
      ilimit=dint(dble(n)*dsqrt(ov/omega)+0.5d0)
      write (6,*) 'ilimit =',ilimit
      write (6,*) 'Ov =',ov
      dr=rmax/(ksiih*dble(n-1))
      df=2d0*pii/dble(m-2)
      d0=dpera/ksiih
      l0=lp
      omega0=omega*ksiih/vd
      ov0=ov*ksiih/vd
c     * LUETAAN ALKUARVAUS *
      if (query.eq.1) then
         do 16 i=1,ilimit
            do 15 j=2,m-1
               fii=(dble(j)-1.5d0)*df
c               fii=(1d0-dexp(-(dble(j-2)/30d0)**2))*pii
c               nv(1,i,j)=-cos(fii)
c               nv(2,i,j)=sin(fii)
c               nv(3,i,j)=0d0
               nv(1,i,j)=0d0
               nv(2,i,j)=0d0
               nv(3,i,j)=1d0
 15         continue
 16      continue
         do 20 i=ilimit+1,n
            nv(1,i,2)=0d0
            nv(2,i,2)=1d0
            nv(3,i,2)=0d0
            nv(1,i,m-1)=0d0
            nv(2,i,m-1)=-1d0
            nv(3,i,m-1)=0d0
            do 19 j=3,m-2
c               fii=(1d0-dexp(-(dble(j-2)/30d0)**2))*pii
               fii=(dble(j)-2d0)*df
c               nv(1,i,j)=-cos(fii)
c               nv(2,i,j)=sin(fii)
c               nv(3,i,j)=0d0
               nv(1,i,j)=1d0
               nv(2,i,j)=0d0
               nv(3,i,j)=0d0
 19         continue
 20      continue
         nv(1,ilimit,2)=-dsqrt(3d0)/2d0
         nv(2,ilimit,2)=0.5d0
         nv(3,ilimit,2)=0d0
         nv(1,ilimit,m-1)=-dsqrt(3d0)/2d0
         nv(2,ilimit,m-1)=-0.5d0
         nv(3,ilimit,m-1)=0d0
         do 21 i=1,n
            nv(1,i,1)=nv(1,i,m-1)
            nv(1,i,m)=nv(1,i,2)
            nv(2,i,1)=nv(2,i,m-1)
            nv(2,i,m)=nv(2,i,2)
            nv(3,i,1)=nv(3,i,m-1)
            nv(3,i,m)=nv(3,i,2)
 21      continue
      else if (query.eq.2) then
         do 26 i=1,n
            do 25 j=1,m 
               read (31,*) nv(1,i,j)
               read (32,*) nv(2,i,j)
               read (33,*) nv(3,i,j)
 25         continue
 26      continue
         nv(1,ilimit,2)=-dsqrt(3d0)/2d0
         nv(2,ilimit,2)=0.5d0
         nv(3,ilimit,2)=0d0
         nv(1,ilimit,m-1)=-dsqrt(3d0)/2d0
         nv(2,ilimit,m-1)=-0.5d0
         nv(3,ilimit,m-1)=0d0
         do 27 i=ilimit+1,n
            nv(1,i,2)=0d0
            nv(2,i,2)=1d0
            nv(3,i,2)=0d0
            nv(1,i,m-1)=0d0
            nv(2,i,m-1)=-1d0
            nv(3,i,m-1)=0d0
 27      continue
         do 28 i=ilimit,n
            nv(1,i,1)=nv(1,i,m-1)
            nv(1,i,m)=nv(1,i,2)
            nv(2,i,1)=nv(2,i,m-1)
            nv(2,i,m)=nv(2,i,2)
            nv(3,i,1)=nv(3,i,m-1)
            nv(3,i,m)=nv(3,i,2)
 28      continue
      else 
         goto 5000
      end if
c      goto 999
      rewind 31
      rewind 32
      rewind 33
      rewind 35
      iter=0
c      goto 1500
c     KOMMENTOI SEURAAVA RIVI POIS KUN ET TARKISTA DERIVOINTEJA
c      goto 3003
      sum=0d0
 100  iter=iter+1
      if (iter.gt.iterm) goto 999
      call matrix(n,m,nv,a)
      write (6,*) iter,energy(n,m,a),dsqrt(sum)
c     PÄÄITERAATIOSILMUKKA:
c     KIERRETÄÄN n:ÄÄ VÄÄNTÖMOMENTIN t MUKAISESTI
      do 250 it=1,itm
         call matrix(n,m,nv,a)
         call gradient(n,m,a,gr)
         call torque(n,m,nv,gr,t)
         sum=0d0
         do 201 i=2,n
            do 200 j=2,m-1
               dnx=(t(2,i,j)*nv(3,i,j)-t(3,i,j)*nv(2,i,j))
               dny=(t(3,i,j)*nv(1,i,j)-t(1,i,j)*nv(3,i,j))
               dnz=(t(1,i,j)*nv(2,i,j)-t(2,i,j)*nv(1,i,j))
               sum=sum+dnx**2+dny**2+dnz**2
               nvf(1,i,j)=nv(1,i,j)-eps*dnx
               nvf(2,i,j)=nv(2,i,j)-eps*dny
               nvf(3,i,j)=nv(3,i,j)-eps*dnz
               nlen=dsqrt(nvf(1,i,j)**2+nvf(2,i,j)**2+nvf(3,i,j)**2)
               nv(1,i,j)=nvf(1,i,j)/nlen
               nv(2,i,j)=nvf(2,i,j)/nlen
               nv(3,i,j)=nvf(3,i,j)/nlen
 200        continue
 201     continue
         dnx=(t(2,1,1)*nv(3,1,1)-t(3,1,1)*nv(2,1,1))
         dny=(t(3,1,1)*nv(1,1,1)-t(1,1,1)*nv(3,1,1))
         dnz=(t(1,1,1)*nv(2,1,1)-t(2,1,1)*nv(1,1,1))
         sum=sum+dnx**2+dny**2+dnz**2
         nvf(1,1,1)=nv(1,1,1)-eps*dnx
         nvf(2,1,1)=nv(2,1,1)-eps*dny
         nvf(3,1,1)=nv(3,1,1)-eps*dnz
         nlen=dsqrt(nvf(1,1,1)**2+nvf(2,1,1)**2+nvf(3,1,1)**2)
         nv(1,1,1)=nvf(1,1,1)/nlen
         nv(2,1,1)=nvf(2,1,1)/nlen
         nv(3,1,1)=nvf(3,1,1)/nlen
         do 221 i=2,n
            nv(1,i,1)=nv(1,i,m-1)
            nv(1,i,m)=nv(1,i,2)
            nv(2,i,1)=nv(2,i,m-1)
            nv(2,i,m)=nv(2,i,2)
            nv(3,i,1)=nv(3,i,m-1)
            nv(3,i,m)=nv(3,i,2)
 221      continue
 250  continue
      goto 100
 999  call matrix(n,m,nv,a)
      do 1001 i=1,n
         do 1000 j=1,m
            write (31,*) nv(1,i,j)
            write (32,*) nv(2,i,j)
            write (33,*) nv(3,i,j)
 1000    continue
 1001 continue
 1500 big=3.5d0*ksiih*10d0
      do 2001 i=2,n,2
c      do 2001 i=2,n
         napu=2
         if (i.le.ilimit) napu=4
         do 2000 j=2,m-1,napu
c         do 2000 j=2,m-1
            fii=(dble(j)-1.5d0)*df
c            xx=dble(i-1)*dr*cos(fii)/rmax*ksiih
c            yy=dble(i-1)*dr*sin(fii)/rmax*ksiih
            xx=dble(i-1)*dr*cos(fii)
            yy=dble(i-1)*dr*sin(fii)
            write (40,*) xx,yy,1d0-nv(3,i,j)**2
            write (35,*) big*xx,big*yy,
     $           cos(fii)*nv(1,i,j)-sin(fii)*nv(2,i,j),
     $           sin(fii)*nv(1,i,j)+cos(fii)*nv(2,i,j)
 2000    continue
 2001 continue
      write (35,*) 0,0,nv(1,1,1),nv(2,1,1)
      call matrix(n,m,nv,a)
      call energy2(n,m,a)
c     TALLENNETAAN TEKSTUURISPEKTRI
      call divide2(n,m,nv,nmr,summa)
      do 2500 i=1,number+1
         nmr(i)=dble(number)*nmr(i)/summa
         write (51,*) dble(i-1)/dble(number),nmr(i)
 2500 continue
c     TALLENNETAAN TEKSTUURISPEKTRI (OSITETTU)
      call kdivide(n,m,nv,knmr,ksumma)
      do 3001 j=1,m
         do 3000 i=1,number+1
            knmr(j,i)=dble(number+1)*knmr(j,i)/ksumma(j)
 3000    continue
 3001 continue
      do 3002 i=1,number+1
         write (56,*) dble(i-1)/dble(number),knmr(int(n/2),i)
 3002 continue
c     * DERIVOINNIN TESTAUSOSA *
c 3003 eps=1d-3
c      call matrix(n,m,nv,a)
c      call gradient(n,m,a,gr)
cc      write (6,*) 'Gradient=',gr(3,3,1,1)
c      do 3100 i=2,n
c         do 3090 j=2,m-1
c            a(2,1,i,j)=a(2,1,i,j)+eps
c            e2=energy(n,m,a)
c            a(2,1,i,j)=a(2,1,i,j)-2d0*eps
c            e1=energy(n,m,a)
c            a(2,1,i,j)=a(2,1,i,j)+eps
c            write (20,*) i,j,gr(2,1,i,j),(e2-e1)/(2d0*eps),
c     $           dabs(gr(2,1,i,j)-(e2-e1)/(2d0*eps))
ccc            write (20,*) i,j,gr(1,1,i,j)
c 3090    continue
c 3100 continue
c      goto 5000
c     TALLENNETAAN LEVENNETTY SPEKTRI
      call convol(nmr,con,vaim1,vaim2)
      do 3500 i=-2*number,2*number
         write (50,*) dble(i)/dble(number),con(i)
 3500 continue
c     TALLENNETAAN LEVENNETTY SPEKTRI (OSITETTU)
      call kconvol(m,knmr,kcon,vaim1,vaim2)
      do 4001 j=1,41
         iso=0d0
         do 4000 i=-2*number,2*number
            write (57,*) dble(i)/dble(number),kcon(5*(j-1),i)+2d0*
     $           dble(j)
            if (kcon(5*(j-1),i).gt.iso) iso=kcon(5*(j-1),i)
 4000    continue
c         write (55,*) (len*5d0)*dble(j-1-20)/dble(20),iso
 4001 continue
      do 4500 i=1,n
c         write (60,*) 1d0-nv(3,i,int(n/2))**2,1d2*
c     $        (dacos(nv(3,i+1,int(n/2)))-dacos(nv(3,i,int(n/2))))/dr
         write (60,*) dble(i-1)/dble(n),
     $        dacos(nv(3,i,int(m/2)))*180d0/pii
 4500 continue
 5000 end

      subroutine matrix(n,m,nv,a)
c     MUODOSTETAAN n-JAKAUMAA VASTAAVA A-MATRIISI
      implicit none 
      integer i,n,nmax,m,j
      parameter (nmax=200)
      double precision a(3,3,nmax,nmax),nv(3,nmax,nmax),e2,e1,energy
      double precision kos,si,f,axx,axy,axz,ayx,ayy,ayz,azx,azy,azz
      kos=-0.25d0
      si=0.25d0*dsqrt(15d0)
      do 40 i=2,n
         do 30 j=1,m
            a(1,1,i,j)=kos+(1d0-kos)*nv(1,i,j)*nv(1,i,j)
            a(2,1,i,j)=(1d0-kos)*nv(2,i,j)*nv(1,i,j)+si*nv(3,i,j)
            a(3,1,i,j)=(1d0-kos)*nv(3,i,j)*nv(1,i,j)-si*nv(2,i,j)
            a(1,2,i,j)=(1d0-kos)*nv(1,i,j)*nv(2,i,j)-si*nv(3,i,j)
            a(2,2,i,j)=kos+(1d0-kos)*nv(2,i,j)*nv(2,i,j)
            a(3,2,i,j)=(1d0-kos)*nv(3,i,j)*nv(2,i,j)+si*nv(1,i,j)
            a(1,3,i,j)=(1d0-kos)*nv(1,i,j)*nv(3,i,j)+si*nv(2,i,j)
            a(2,3,i,j)=(1d0-kos)*nv(2,i,j)*nv(3,i,j)-si*nv(1,i,j)
            a(3,3,i,j)=kos+(1d0-kos)*nv(3,i,j)*nv(3,i,j)
 30      continue
 40   continue
      axx=kos+(1d0-kos)*nv(1,1,1)*nv(1,1,1)
      ayx=(1d0-kos)*nv(2,1,1)*nv(1,1,1)+si*nv(3,1,1)
      azx=(1d0-kos)*nv(3,1,1)*nv(1,1,1)-si*nv(2,1,1)
      axy=(1d0-kos)*nv(1,1,1)*nv(2,1,1)-si*nv(3,1,1)
      ayy=kos+(1d0-kos)*nv(2,1,1)*nv(2,1,1)
      azy=(1d0-kos)*nv(3,1,1)*nv(2,1,1)+si*nv(1,1,1)
      axz=(1d0-kos)*nv(1,1,1)*nv(3,1,1)+si*nv(2,1,1)
      ayz=(1d0-kos)*nv(2,1,1)*nv(3,1,1)-si*nv(1,1,1)
      azz=kos+(1d0-kos)*nv(3,1,1)*nv(3,1,1)
c      azz=azz+1d-3
      do 100 j=1,m
         f=(dble(j)-1.5d0)*2d0*3.14159265d0/dble(m-2)
         a(1,1,1,j)=axx*cos(f)**2+sin(f)*cos(f)*(axy+ayx)+ayy*sin(f)**2
         a(2,1,1,j)=ayx*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-axy*sin(f)**2
         a(3,1,1,j)=azx*cos(f)+azy*sin(f)
         a(1,2,1,j)=axy*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-ayx*sin(f)**2
         a(2,2,1,j)=ayy*cos(f)**2-sin(f)*cos(f)*(axy+ayx)+axx*sin(f)**2
         a(3,2,1,j)=azy*cos(f)-azx*sin(f)
         a(1,3,1,j)=axz*cos(f)+ayz*sin(f)
         a(2,3,1,j)=ayz*cos(f)-axz*sin(f)
         a(3,3,1,j)=azz
 100  continue
c      e2=energy(n,m,a)
c      azz=azz-2d-3
c      do 200 j=1,m
c         f=(dble(j)-1.5d0)*2d0*3.14159265d0/dble(m-2)
c         a(1,1,1,j)=axx*cos(f)**2+sin(f)*cos(f)*(axy+ayx)+ayy*sin(f)**2
c         a(2,1,1,j)=ayx*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-axy*sin(f)**2
c         a(3,1,1,j)=azx*cos(f)+azy*sin(f)
c         a(1,2,1,j)=axy*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-ayx*sin(f)**2
c         a(2,2,1,j)=ayy*cos(f)**2-sin(f)*cos(f)*(axy+ayx)+axx*sin(f)**2
c         a(3,2,1,j)=azy*cos(f)-azx*sin(f)
c         a(1,3,1,j)=axz*cos(f)+ayz*sin(f)
c         a(2,3,1,j)=ayz*cos(f)-axz*sin(f)
c         a(3,3,1,j)=azz
c 200  continue
c      e1=energy(n,m,a)
c      write (6,*) 'Quotient=',(e2-e1)/2d-3
c      azz=azz+1d-3
c      do 300 j=1,m
c         f=(dble(j)-1.5d0)*2d0*3.14159265d0/dble(m-2)
c         a(1,1,1,j)=axx*cos(f)**2+sin(f)*cos(f)*(axy+ayx)+ayy*sin(f)**2
c         a(2,1,1,j)=ayx*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-axy*sin(f)**2
c         a(3,1,1,j)=azx*cos(f)+azy*sin(f)
c         a(1,2,1,j)=axy*cos(f)**2-sin(f)*cos(f)*(axx-ayy)-ayx*sin(f)**2
c         a(2,2,1,j)=ayy*cos(f)**2-sin(f)*cos(f)*(axy+ayx)+axx*sin(f)**2
c         a(3,2,1,j)=azy*cos(f)-azx*sin(f)
c         a(1,3,1,j)=axz*cos(f)+ayz*sin(f)
c         a(2,3,1,j)=ayz*cos(f)-axz*sin(f)
c         a(3,3,1,j)=azz
c 300  continue
      return
      end

      subroutine torque(n,m,nv,gr,t)
c     LASKETAAN E-L:N YHTÄLÖISTÄ SEURAAVA n:N VÄÄNTÖMOMENTTI
      implicit none 
      integer i,n,nmax,m,j
      parameter (nmax=200)
      double precision nv(3,nmax,nmax),gr(3,3,nmax,nmax)
      double precision t(3,nmax,nmax)
      double precision kos,si
      kos=-0.25d0
      si=0.25d0*dsqrt(15d0)
      do 300 i=2,n
         do 200 j=2,m-1
            t(1,i,j)=nv(2,i,j)*(1d0-kos)*(nv(1,i,j)*(gr(1,3,i,j)+
     $           gr(3,1,i,j))+nv(2,i,j)*(gr(2,3,i,j)+gr(3,2,i,j))+
     $           2d0*nv(3,i,j)*gr(3,3,i,j))-nv(2,i,j)*si*(gr(1,2,i,j)
     $           -gr(2,1,i,j))-nv(3,i,j)*(1d0-kos)*
     $           (nv(1,i,j)*(gr(1,2,i,j)+gr(2,1,i,j))
     $           +nv(3,i,j)*(gr(2,3,i,j)+gr(3,2,i,j))+
     $           2d0*nv(2,i,j)*gr(2,2,i,j))
     $           +nv(3,i,j)*si*(gr(3,1,i,j)-gr(1,3,i,j))
            t(2,i,j)=nv(3,i,j)*(1d0-kos)*(nv(2,i,j)*(gr(1,2,i,j)+
     $           gr(2,1,i,j))+nv(3,i,j)*(gr(1,3,i,j)+gr(3,1,i,j))+
     $           2d0*nv(1,i,j)*gr(1,1,i,j))-nv(3,i,j)*si*(gr(2,3,i,j)
     $           -gr(3,2,i,j))-nv(1,i,j)*(1d0-kos)*(nv(1,i,j)*
     $           (gr(1,3,i,j)+gr(3,1,i,j))+nv(2,i,j)*(gr(2,3,i,j)+
     $           gr(3,2,i,j))+2d0*nv(3,i,j)*gr(3,3,i,j))
     $           +nv(1,i,j)*si*(gr(1,2,i,j)-gr(2,1,i,j))
            t(3,i,j)=nv(1,i,j)*(1d0-kos)*(nv(1,i,j)*(gr(1,2,i,j)+
     $           gr(2,1,i,j))+nv(3,i,j)*(gr(2,3,i,j)+gr(3,2,i,j))+
     $           2d0*nv(2,i,j)*gr(2,2,i,j))-nv(1,i,j)*si*(gr(3,1,i,j)
     $           -gr(1,3,i,j))-nv(2,i,j)*(1d0-kos)*(nv(2,i,j)*
     $           (gr(2,1,i,j)+gr(1,2,i,j))+nv(3,i,j)*(gr(3,1,i,j)+
     $           gr(1,3,i,j))+2d0*nv(1,i,j)*gr(1,1,i,j))
     $           +nv(2,i,j)*si*(gr(2,3,i,j)-gr(3,2,i,j))
 200     continue
 300  continue
      t(1,1,1)=nv(2,1,1)*(1d0-kos)*(nv(1,1,1)*(gr(1,3,1,1)+
     $     gr(3,1,1,1))+nv(2,1,1)*(gr(2,3,1,1)+gr(3,2,1,1))+
     $     2d0*nv(3,1,1)*gr(3,3,1,1))-nv(2,1,1)*si*(gr(1,2,1,1)
     $     -gr(2,1,1,1))-nv(3,1,1)*(1d0-kos)*
     $     (nv(1,1,1)*(gr(1,2,1,1)+gr(2,1,1,1))
     $     +nv(3,1,1)*(gr(2,3,1,1)+gr(3,2,1,1))+
     $     2d0*nv(2,1,1)*gr(2,2,1,1))
     $     +nv(3,1,1)*si*(gr(3,1,1,1)-gr(1,3,1,1))
      t(2,1,1)=nv(3,1,1)*(1d0-kos)*(nv(2,1,1)*(gr(1,2,1,1)+
     $     gr(2,1,1,1))+nv(3,1,1)*(gr(1,3,1,1)+gr(3,1,1,1))+
     $     2d0*nv(1,1,1)*gr(1,1,1,1))-nv(3,1,1)*si*(gr(2,3,1,1)
     $     -gr(3,2,1,1))-nv(1,1,1)*(1d0-kos)*(nv(1,1,1)*
     $     (gr(1,3,1,1)+gr(3,1,1,1))+nv(2,1,1)*(gr(2,3,1,1)+
     $     gr(3,2,1,1))+2d0*nv(3,1,1)*gr(3,3,1,1))
     $     +nv(1,1,1)*si*(gr(1,2,1,1)-gr(2,1,1,1))
      t(3,1,1)=nv(1,1,1)*(1d0-kos)*(nv(1,1,1)*(gr(1,2,1,1)+
     $     gr(2,1,1,1))+nv(3,1,1)*(gr(2,3,1,1)+gr(3,2,1,1))+
     $     2d0*nv(2,1,1)*gr(2,2,1,1))-nv(1,1,1)*si*(gr(3,1,1,1)
     $     -gr(1,3,1,1))-nv(2,1,1)*(1d0-kos)*(nv(2,1,1)*
     $     (gr(2,1,1,1)+gr(1,2,1,1))+nv(3,1,1)*(gr(3,1,1,1)+
     $     gr(1,3,1,1))+2d0*nv(1,1,1)*gr(1,1,1,1))
     $     +nv(2,1,1)*si*(gr(2,3,1,1)-gr(3,2,1,1))
      return
      end

      subroutine divide(n,m,nv,nmr,summa)
c     MUODOSTETAAN sin^2 beta -JAKAUMA
      implicit none
      integer i,j,n,nmax,luku,number,m
      parameter (nmax=200,number=100)
      double precision nmr(number+1)
      double precision nv(3,nmax,nmax),apu,summa
      summa=0d0
      do 10 i=1,number+1
         nmr(i)=0d0
 10   continue
      do 100 i=1,n
         do 99 j=1,m
            apu=dble(number+1)*(1d0-nv(3,i,j)**2)
            luku=dint(apu)
            if (luku.gt.number) luku=number
            nmr(luku+1)=nmr(luku+1)+dble(i-1)
 99      continue
 100  continue
      do 150 i=1,number+1
         summa=summa+nmr(i)
 150  continue
      return 
 200  end

      subroutine divide2(n,m,nv,nmr,summa)
      implicit none
      integer i,j,n,m,nmax,luku,number,niso,k,nk
      parameter (nmax=200,number=100,niso=40000)
      double precision nmr(number+1),maksi,kerr
      double precision nv(3,nmax,nmax),apu,summa,uusi(niso,nmax)
      summa=0d0
      maksi=0d0
      nk=500
      do 10 i=1,number+1
         nmr(i)=0d0
 10   continue
c     *******     uusi = cos(beta)    ******************
c     *******     100*tihennettynä     ******************
      do 30 i=1,n-1
         do 29 j=2,m-1
            uusi(nk*(i-1)+1,j)=nv(3,i,j)
            do 20 k=1,nk-1
               uusi(nk*(i-1)+1+k,j)=nv(3,i,j)+1d0/dble(nk)*
     $              dble(k)*(nv(3,i+1,j)-nv(3,i,j))
 20         continue
 29      continue
 30   continue
      do 40 j=2,m-1
         uusi(nk*n+1,j)=nv(3,n+1,j)
 40   continue
c     ******* Jaetaan laatikoihin  *****************
      do 100 i=1,nk*(n-1)+1
         do 99 j=2,m-1
            apu=dble(number+1)*(1d0-uusi(i,j)**2)
            luku=dint(apu)
            if (luku.gt.number) luku=number
c            nmr(luku+1)=nmr(luku+1)+dble(i-1)/dble(nk*(n-1))
            nmr(luku+1)=nmr(luku+1)+dble(i-1)
 99      continue
 100  continue
      do 150 i=1,number+1
         kerr=1d0
         if ((i.eq.1).or.(i.eq.(number+1))) kerr=0.5d0
         summa=summa+kerr*nmr(i)
 150  continue
      return 
 200  end

      subroutine convol(nmr,con,vaim1,vaim2)
c     LASKETAAN LEVENEMÄKONVOLUUTIO
      implicit none
      integer i,j,number,nmax,sim
      double precision pii
      parameter (number=100,pii=3.14159265359d0)
      double precision nmr(number+1),gamma,apu,vaim1,kerr
      double precision con(-2*number:2*number),delx,dely,vaim2
      delx=1d0/dble(number)
      dely=delx
      do 10 j=-2*number,2*number
         con(j)=0d0
         do 5 i=1,number+1
c     *********  Trapetsi / Simpson **********************
            gamma=(vaim1+vaim2*dble(i-1)/dble(number))/2d0
            apu=gamma*dely/pii
c            apu=gamma*dely/(3d0*pii)
c            if ((i.eq.1).or.(i.eq.number)) then 
c               kerr=1d0
c               sim=2
c            else
c               if (sim.eq.2) then
c                  kerr=4d0
c                  sim=4
c               else
c                  kerr=2d0
c                  sim=2
c               end if
c            end if            
            if ((i.eq.1).or.(i.eq.number+1)) then 
               kerr=0.5d0
            else
               kerr=1d0
            end if
            con(j)=con(j)+kerr*apu*nmr(i)/((delx*dble(j)-
     $           dely*dble(i-1))**2+gamma**2)
 5       continue
 10   continue
      return
      end

      subroutine sconvol(nmr,con,vaim1,vaim2)
c     LASKETAAN (SPIND.) LEVENEMÄKONVOLUUTIO
      implicit none
      integer i,j,number,nmax,sim
      double precision pii,nyyb,nyyl,dkk,dc
      parameter (number=100,pii=3.14159265359d0)
      double precision nmr(number+1),gamma,apu,vaim1,kerr,dif,sb
      double precision con(-2*number:2*number),delx,dely,vaim2
      delx=1d0/dble(number)
      dely=delx
      nyyb=60d3
      nyyl=688d3
      dkk=3d-5*(1.1d0/2d-3)**2
      dc=2d0*dkk/(pii*nyyb**2/nyyl)
      write (6,*) 'Diffusion constant =',dc
      do 10 j=-2*number,2*number
         con(j)=0d0
         do 5 i=1,number+1
c     *********  Trapetsi / Simpson **********************
            sb=dble(i-1)*dely
            if ((dble(i-1)/dble(number)).lt.(0.77d0)) then
               dc=2d0*dkk/(pii*nyyb**2/nyyl)
            else
               dc=0d0
            end if
            gamma=(vaim1+vaim2*dble(i-1)/dble(number))
            apu=dely/(3d0*pii)
            if ((i.eq.1).or.(i.eq.number)) then 
               kerr=1d0
               sim=2
            else
               if (sim.eq.2) then
                  kerr=4d0
                  sim=4
               else
                  kerr=2d0
                  sim=2
               end if
            end if
            dif=delx*dble(j)-dely*dble(i-1)
            con(j)=con(j)+kerr*apu*nmr(i)*(gamma/(dif**2+gamma**2)
     $           -dc*(4d0*sb*(1d0-sb)+dif*(1d0-2d0*sb))*
     $           (dif**4-6d0*(dif*gamma)**2+gamma**4)/
     $           (dif**2+gamma**2)**4)
 5       continue
 10   continue
      return
      end

      subroutine kdivide(n,m,nv,knmr,ksumma)
c     MUODOSTETAAN sin^2 beta -JAKAUMA
      implicit none
      integer i,j,n,nmax,luku,number,m
      parameter (nmax=200,number=100)
      double precision knmr(m,number+1)
      double precision nv(3,nmax,nmax),apu,ksumma(m)
      do 11 j=1,m
         ksumma(j)=0d0
         do 10 i=1,number+1
            knmr(j,i)=0d0
 10      continue
 11   continue
      do 100 i=1,m
         do 99 j=1,n
            apu=dble(number+1)*(1d0-nv(3,j,i)**2)
            luku=dint(apu)
            if (luku.gt.number) luku=number
            knmr(i,luku+1)=knmr(i,luku+1)+dble(j-1)
 99      continue
 100  continue
      do 151 j=1,m
         do 150 i=1,number+1
            ksumma(j)=ksumma(j)+knmr(j,i)
 150     continue
 151  continue
      return 
 200  end

      subroutine kconvol(m,knmr,kcon,vaim1,vaim2)
c     LASKETAAN LEVENEMÄKONVOLUUTIO
      implicit none
      integer i,j,k,number,nmax,sim,m
      double precision pii
      parameter (number=100,pii=3.14159265359d0)
      double precision knmr(m,number+1),gamma,apu,vaim1,kerr
      double precision kcon(m,-2*number:2*number),delx,dely,vaim2
      delx=1d0/dble(number)
      dely=delx
      do 20 k=1,m
         do 10 j=-2*number,2*number
            kcon(k,j)=0d0
            do 5 i=1,number+1
c     *********  Trapetsi / Simpson **********************
               gamma=(vaim1+vaim2*dble(i-1)/dble(number))/2d0
c     apu=gamma*dely/(3d0*pii)
               apu=gamma*dely/pii
c     if ((i.eq.1).or.(i.eq.number+1)) then 
c     kerr=1d0
c     sim=2
c     else
c     if (sim.eq.2) then
c     kerr=4d0
c     sim=4
c     else
c     kerr=2d0
c     sim=2
c     end if
c     end if            
               if ((i.eq.1).or.(i.eq.number+1)) then 
                  kerr=0.5d0
               else
                  kerr=1d0
               end if
               kcon(k,j)=kcon(k,j)+kerr*apu*knmr(k,i)/((delx*dble(j)-
     $              dely*dble(i-1))**2+gamma**2)
 5          continue
 10      continue
 20   continue
      return
      end
