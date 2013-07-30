      program area
      implicit none 
      integer i,j,m,n,number
      double precision dx
      parameter (number=200,dx=1d-2)
      double precision x,y,sum
      open (18,file='testi')
      sum=0d0
      do 100 i=-number,number
         read (18,*) x,y
         if ((y.gt.0d0).and.(x.gt.0.5d0)) sum=sum+y*dx
 100  continue
      write (6,*) sum
      end

      
