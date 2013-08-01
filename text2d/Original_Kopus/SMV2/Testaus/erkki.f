c main program for testing
      subroutine vakiot(temp,omega,p,h,rtot,dpera,vd,lp,ksiih,de,
     $     dhh,tau,nyyb,vaim1,vaim2)
      implicit none
      integer i,nfreq,iparam
      double precision temp,gap,gasc,pii,ap,ghv,mbare,rho,h,vaim2
      double precision p,f1s,f0s,f0a,f1a,vol,tc,sheat,dcpcb,dcpca
      double precision gd,kb,hbar,gamm,khi,y,g3,g5,g7,vd,tau
      double precision de,gg,rs,dp,xi,kf,omega,lp,rvor,ksiih,rtot
      double precision ffunc,vaim1,f,lam,khikhi0,dpera
      double precision nyyb,nyyl,dens,vaim3,dhh
      parameter (pii=3.1415926535898d0,kb=1.38d-16,hbar=1.05d-27)
      parameter (gamm=-2.04d4,mbare=5.01d-24,gasc=8.31d7)
      iparam=0
      nyyl=688d3*h/200d0
c*****    Magneettikent‰n ep‰homogeenisuus  ********************* 
      vaim1=dhh*4d0*(nyyl/nyyb)**2
c      write (6,*) 'DeltaW / WL=',0.5d0*(nyyb/nyyl)**2*vaim1
      rvor=dsqrt(hbar/(2d0*mbare*omega))
      if (omega.eq.0d0) rvor=0d0
      call ferpar(p,f1s,f0s,f0a,f1a,vol,tc,sheat,dcpcb,dcpca,
     $   gd,iparam)
c      dens=3d0*gasc*sheat/(2d0*pii*pii*kb*kb*vol)
      rho=(mbare/vol)*6.02d23
      kf=(3d0*pii**2*6.02d23/vol)**(1d0/3d0)
      dens=kf*mbare*(1d0+f1s/3d0)/(2d0*pii**2*hbar**2)
      call bcsgap(temp,gap)
      call wcp(temp,dcpcb,gap)
      write (6,*) 'Gap/Tc =',gap
      g3=0d0
      g5=0d0
      g7=0d0
      do 50 i=-200,200
         g3=g3+pii*gap**2*temp*((pii*temp*(2d0*dble(i)+1d0))**2
     $        +gap*gap)**(-1.5d0)
         g5=g5+pii*gap**4*temp*((pii*temp*(2d0*dble(i)+1d0))**2
     $        +gap*gap)**(-2.5d0)
         g7=g7+pii*gap**6*temp*((pii*temp*(2d0*dble(i)+1d0))**2
     $        +gap*gap)**(-3.5d0)
 50   continue
      y=1d0-g3
c*****      Leggett-Takagi -relaksaatio    ********************* 
      f=ffunc(temp,gap)
      lam=2d0*(1d0-f)/(2d0+y)
      khikhi0=1d0/(1d0+f0a*(2d0/3d0+y/3d0))
      vaim2=(1d0-lam)/lam*khikhi0*tau*4d0*pii*nyyb**2/nyyl
      vaim3=(1d0-lam)/lam*khikhi0*tau*4d0*pii**2*nyyb**4/nyyl**2
c      write (6,*) 'Kappa =',tau*khikhi0*(1d0-lam)/lam
      write (6,*) 'LT-Vaim =',vaim3
      write (6,*) 'dH-Vaim =',vaim1*pii*nyyb**2/nyyl
c*****     Vaimennusosuus loppuu          *********************
      khi=-2d0*dens*(hbar*gamm/2d0)**2*((2d0/3d0+y/3d0)/(1d0+
     $     f0a*(2d0+y)/3d0)-1d0/(1d0+f0a))
      ap=2.5d32*gd*(0.5d0*hbar*gamm/(1d0+f0a*(2d0+y)/3d0))**2*
     $     (5d0-3d0*g5/g3)
      ghv=(rho*(1d0+f1s/3d0)/(gap*kb*tc*1d-3*
     $     (1d0+f1s*y/3d0))**2)*
     $     (0.5d0*hbar*gamm/(1d0+f0a*(2d0+y)/3d0))**2
     $     *(g3-0.9d0*g5+0.9d0*g5**2/g3-1.5d0*g7)
      gg=hbar**2*rho/(4d1*mbare**2*(1d0+f1s/3d0))*(1d0+f1a/3d0)
     $     *(1d0-y)/(1d0+f1a*(2d0+3d0*y)/15d0)
      de=f1a/3d0*(1d0-y)/(1d0+f1a*y/3d0)
      xi=hbar**2*kf/(mbare*kb*tc*1d-3*dsqrt(1d1)
     $     *gap*(1d0+f1s/3d0))
      rs=2.5d0*xi
      dp=0.5d0*khi*xi*2.2d0
      lp=5d0*ghv*hbar*omega*(dlog(rvor/xi)-0.75d0)/(2d0*mbare*ap)
      ksiih=dsqrt(65d0*gg/(8d0*ap*h**2))
      vd=dsqrt(2d0*ap/(5d0*ghv))
      dpera=dp/ap
      return
      end

c Subroutine to calculate the bcs gap function in the pure case (SC or 3He-B)
c Essentially an exact calculation, see (Ajop6,10)
c The idea is to sum exactly the first, say 10, terms in the matsubara
c sum, and estimate the rest by an integral and a Euler-McLaurin 
c correction term
c input:
c   temp= temperature/T_c
c output:
c   gap= the energy gap/T_c
      subroutine bcsgap(temp,gap)
      implicit none
      double precision temp,gap,pii2,tn
      double precision fprime,eps,sq,gaps,f,x
      integer iter,ifreq,nsplit
      parameter(pii2=2d0*3.1415926535897932d0,nsplit=30)
      iter=0
      if(temp.ge.1d0)then
         gap=0d0
	 goto 200
      end if
      gap=1.7638d0*dsqrt(1d0-temp)/pii2
 10   iter=iter+1
      gaps=gap*gap
      tn=temp*dble(nsplit)
      sq=dsqrt(tn**2+gaps)
      f=dlog((tn+sq)/(2d0*dble(nsplit)))
     $   -(1d0/dble(nsplit)**2-dble(nsplit)*(temp/sq)**3)/24d0
      fprime=gap/(sq*(tn+sq))-tn*temp**2*gap/(8d0*sq**5)
      do 100 ifreq=1,nsplit
         eps=dble(ifreq)-0.5d0
	 sq=dsqrt((temp*eps)**2+gaps)
	 f=f+1d0/eps-temp/sq
	 fprime=fprime+temp*gap/sq**3
 100  continue
      x=f/fprime
      gap=gap-x
      if(abs(x).ge.1d-8/pii2) goto 10
      gap=pii2*gap
 200  continue
C      write(6,1)gap,iter,temp
 1    format(' gap=',f8.6,' iter=',i3,' for temp=',f8.6)
      return
      end


c laskee matsubara cuf-off energian kun on annettu tarkkuus
c accur, jolla PUHTAAN nesteen gap pit‰‰ olla on oikein.
c input:
c   temp= temperature/Tc
c   accur= allowed error in the gap/Tc
c output:
c   nfreq= number of positive matsubara energies needed to achieve
c          the desired accuracy (assuming that all higher energies are
c          simply omitted from the sum in the gap equation)
      subroutine cutoff(temp,accur,nfreq)
      implicit none
      integer nfreq
      double precision temp, gap, fprime, err
      double precision pi, zeta37, eps, sqterm
      double precision acc, accur,small
      parameter (pi = 3.141592654d0,small=.000001d0)
      parameter (zeta37 = 1.2020569d0*7d0)
      if(accur.ge.0.5d0)then
         nfreq=int(accur+small)
      else
         acc = accur/(2d0*pi*temp)
         gap = 1.7638d0*dsqrt(1d0-temp)/(2d0*pi*temp)
         fprime = 0d0
         err = zeta37
         nfreq = 0
 80      nfreq = nfreq + 1
         eps = dble(nfreq) - 0.5d0
         err = err - 1d0/eps**3
         sqterm = dsqrt(eps**2 + gap**2)
c seuraava 2:nen tulee koska virhe itse asiassa on gap**2*err/2
         fprime = fprime + 2d0/sqterm**3
         if (gap*err/fprime .gt. acc) goto 80
      end if
      return
      end


c subroutine to calculate the trivial strong-coupling corretion
c for B phase (according to the table in Serene-Rainer review)
c uses linear extrapolation for large dcpcn
c input
c   gap = bare energy gap of the B phase (unit=k_BT_c)
c   dcpcn=specific heat jump in B phase =deltaC/Cnormal at T=T_c
c   temp = temperature/superfluid transition temperature
c output:
c   gap= corrected gap 
      subroutine wcp(temp,dcpcn,gap)
      implicit none
      integer j,it,itp,icp,ic
      double precision del,corr
      double precision temp,dcpcn,gap,wt1,wt2,wc1,wc2
      double precision c(5),x(0:10,5)
      data (c(j),j=1,5) /1.43d0,1.6d0,1.8d0,2d0,2.2d0/
      data (x(10,j),j=1,5) /1d0,1.056d0,1.115d0,1.171d0,1.221d0/
      data (x(9,j),j=1,5) /1d0,1.048d0,1.097d0,1.141d0,1.18d0/
      data (x(8,j),j=1,5) /1d0,1.041d0,1.083d0,1.119d0,1.15d0/
      data (x(7,j),j=1,5) /1d0,1.036d0,1.072d0,1.102d0,1.128d0/
      data (x(6,j),j=1,5) /1d0,1.032d0,1.063d0,1.089d0,1.112d0/
      data (x(5,j),j=1,5) /1d0,1.028d0,1.056d0,1.079d0,1.099d0/
      data (x(4,j),j=1,5) /1d0,1.026d0,1.051d0,1.073d0,1.091d0/
      data (x(3,j),j=1,5) /1d0,1.024d0,1.049d0,1.069d0,1.086d0/
      data (x(2,j),j=1,5) /1d0,1.024d0,1.048d0,1.068d0,1.085d0/
      data (x(1,j),j=1,5) /1d0,1.024d0,1.048d0,1.068d0,1.085d0/
      data (x(0,j),j=1,5) /1d0,1.024d0,1.048d0,1.068d0,1.085d0/
      it=int(temp*10d0-0.00001d0)
      if(it.ge.10)return
      itp=it+1
      wt1=(0.1d0*dble(itp)-temp)/0.1d0
      wt2=1d0-wt1
      do 100 ic=1,4
         icp=ic+1
 100  if(dcpcn.lt.c(icp)) goto 110
      ic=4
      icp=5
 110  del=c(icp)-c(ic)
      wc1=(c(icp)-dcpcn)/del
      wc2=1d0-wc1
      corr=wt1*(wc1*x(it,ic)+wc2*x(it,icp))
     $   +wt2*(wc1*x(itp,ic)+wc2*x(itp,icp))
      gap=gap*corr
      return
      end


c This subprogram gives the fermi-liquid parameters as a function of
c pressure.
c input:
c   p = pressure (unit=bar)
c   iparam=0: according to the table of vollhardt and wofle,
c              h.h. hensley, jltp 90, 149 (1993)(exept melting pressure).
c         =1: polynomial fit by Halperin and Varoguaux
c         =2: Fetter values at the melting pressure
c output:
c   f1s,f0s,f0a,f1a = fermi liquid parameters
c   vol= volume per mole (cm^3)
c   tc = transition temperature (mK)
c   sheat=normal state specific heat coefficient C/nRT (unit=1/K)
c   dcpcb=specific heat jump in B phase =deltaC/Cnormal at T=T_c
c   dcpca=specific heat jump in A phase =deltaC/Cnormal at T=T_c
c   gd=dipole coefficient (non T-dependent) (unit=e32 1/erg cm^3)
      subroutine ferpar(p,f1s,f0s,f0a,f1a,vol,tc,sheat,dcpcb,dcpca,
     $   gd,iparam)
      implicit none
      integer j,i,iparam,ip
      double precision p,f1s,f0s,f0a,f1a,vol,tc,sheat,dcpcb
      double precision dcpca,gd,w2,w1,delta
      double precision pt(16),x(16,6)
      data(pt(i),i=1,16)/0d0,1d0,2d0,3d0,6d0,9d0,12d0,15d0,18d0,21d0,
     $   24d0,27d0,30d0,33d0,34.39d0,100000d0/
      data(x(1,j),j=1,6)/ 5.39d0, 9.3d0,-0.698d0,36.84d0,
     $   0.929d0,2.78d0/
      data(x(2,j),j=1,6)/ 5.78d0,11.53d0,-0.710d0,35.74d0,
     $   1.061d0,2.85d0/
      data(x(3,j),j=1,6)/ 6.14d0,13.76d0,-0.718d0,34.78d0,
     $   1.181d0,2.92d0/
      data(x(4,j),j=1,6)/ 6.49d0,15.99d0,-0.724d0,33.95d0,
     $   1.29d0,2.98d0/
      data(x(5,j),j=1,6)/ 7.45d0,22.49d0,-0.734d0,32.03d0,
     $   1.56d0,3.16d0/
      data(x(6,j),j=1,6)/ 8.31d0,29d0,-0.741d0,30.71d0,
     $   1.769d0,3.32d0/
      data(x(7,j),j=1,6)/ 9.09d0,35.42d0,-0.746d0,29.71d0,
     $   1.934d0,3.48d0/
      data(x(8,j),j=1,6)/ 9.85d0,41.73d0,-0.751d0,28.89d0,
     $   2.067d0,3.62d0/
      data(x(9,j),j=1,6)/10.6d0,48.46d0,-0.754d0,28.18d0,
     $   2.177d0,3.77d0/
      data(x(10,j),j=1,6)/11.34d0,55.2d0,-0.756d0,27.55d0,
     $   2.267d0,3.92d0/
      data(x(11,j),j=1,6)/12.07d0,62.16d0,-0.757d0,27.01d0,
     $   2.339d0,4.06d0/
      data(x(12,j),j=1,6)/12.79d0,69.43d0,-0.757d0,26.56d0,
     $   2.395d0,4.21d0/
      data(x(13,j),j=1,6)/13.5d0,77.02d0,-0.754d0,26.17d0,
     $   2.438d0,4.36d0/
      data(x(14,j),j=1,6)/14.21d0,84.79d0,-0.755d0,25.75d0,
     $   2.474d0,4.5d0/
      data(x(15,j),j=1,6)/14.56d0,88.47d0,-0.753d0,25.5d0,
     $   2.491d0,4.56d0/
      data(x(16,j),j=1,6)/14.56d0,88.47d0,-0.753d0,25.5d0,
     $   2.491d0,4.56d0/
      if (p.gt.34.39d0) p=34.39d0
      f1a=-0.5678d0+p*(-0.04753d0+p*(1.791d-3-p*2.273d-5))
      gd=0.27733d0+5.8087d-4*p+2.515d-4*p**2
      if(iparam.eq.0)then
         do 100 i=1,15
            ip=i+1
 100     if(p.lt.pt(ip)) goto 110
 110     delta=pt(ip)-pt(i)
         w1=(pt(ip)-p)/delta
         w2=1d0-w1
         f1s=w1*x(i,1)+w2*x(ip,1)
         f0s=w1*x(i,2)+w2*x(ip,2)
         f0a=w1*x(i,3)+w2*x(ip,3)
         vol=w1*x(i,4)+w2*x(ip,4)
         tc=w1*x(i,5)+w2*x(ip,5)
	 sheat=w1*x(i,6)+w2*x(ip,6)
      else if(iparam.eq.1)then
	 vol=(36.82d0+p*(-1.185d0+p*(0.08609d0+p*(-0.004187d0+p*
     $      (1.071d0-p*1.082d-6)))))
	 tc=0.9294d0+p*(0.1387d0+p*(-0.00693d0+p*(0.0002569d0+
     $      p*(-0.000005725d0+p*0.00000005301d0))))
	 sheat=2.78d0+p*(0.07151d0+p*(-0.001678d0+p*
     $      (0.00005342d0-p*0.000000616d0)))
	 f0s=1000000d0
	 f1s=3d0*(2.8d0+p*(0.1292d0+p*(-0.003188d0+p*
     $      (0.00009372d0-p*0.00000103d0)))-1d0)
	 f0a=-0.7007d0+p*(-0.006232d0+p*(0.0002057d0-
     $      p*0.000001823d0))
      else
c Parmeters used by Fetter in calulation of A-phase hydrostatic
c parameters (tc, p, and f0s not needed, vol and sheat adjusted to give
c the same susceptibility, and gd to give the same rld at temp=0.7)
         p=34.39d0
         vol=25.54d0
         tc=0.9294d0+p*(0.1387d0+p*(-0.00693d0+p*(0.0002569d0+
     $      p*(-0.000005725d0+p*0.00000005301d0))))
         sheat=4.8917914d0
         f0s=1000000d0
         f1s=15.66d0
         f0a=-0.7375d0
         f1a=-0.54d0
         gd=0.85764884d0
      end if
      dcpcb=41.9d0/vol+0.322d0
      dcpca=94.2d0/vol-1.58d0
      return
      end

      double precision function ffunc(temp,gap)
      implicit none
      integer i,ni,niv
      double precision temp,gap,e,isoe,apu1,apu2
      ni=200
      niv=200
      ffunc=0d0
      do 100 i=0,ni*niv
         e=dble(i)/dble(niv)
         isoe=dsqrt(gap**2+e**2)
         apu1=1d0/(2d0*temp)*isoe
         apu2=4d0/(dexp(apu1)+dexp(-apu1))**2
         ffunc=ffunc+1d0/(2d0*temp*dble(niv))*(e/isoe)**2
     $        *apu2
 100  continue
      return
      end
