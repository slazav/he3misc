      program reduc
      implicit none
      integer i,n
      parameter(n=100)
      double precision x1,x2,x,y1,y2,sum1,sum2,sum3
      double precision apu(2*n+1),kerr
      open (15,file='nosmv.dat')
      open (16,file='smv.dat')
      sum1=0d0
      sum2=0d0
c      do 20 i=1,201
c         read(15,*) x,apu(i)
c 20   continue
c      do 50 i=1,2*n
c         x1=-2d0+(i-1)*0.005d0
c         x2=-2d0+i*0.005d0
c         write(18,*) x1,apu(i)
c         write(18,*) x2,0.5d0*(apu(i)+apu(i+1))
c 50   continue
c      write(18,*) 2d0,apu(2*n+1)
c      rewind 18
      do 100 i=1,4*n+1
         kerr=1d0
         if ((i.eq.1).or.(i.eq.(4*n+1))) kerr=0.5d0
         read(15,*) x,y1
         read(16,*) x,y2
         write(17,*) x,y2-y1
         sum1=sum1+kerr*y1
         sum2=sum2+kerr*y2
 100  continue
      sum3=0d0
      rewind 17
      do 200 i=1,4*n+1
         read(17,*) x,y1
         sum3=sum3+y1
 200  continue
      write (6,*) 'nosmv',sum1/dble(n)
      write (6,*) 'smv',sum2/dble(n)
      write (6,*) 'alkuper',sum3/dble(n)
      end

