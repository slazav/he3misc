      program tupla
      implicit none
      integer i,n,m,nmax,j,ilimit,nuusi
      parameter (n=15,m=50,nmax=500,ilimit=8)
      double precision nv(3,nmax,nmax)
      double precision nv2(3,nmax,nmax)
      double precision nv1(3,nmax,nmax),d
      nuusi=2*n-1
      do 100 i=1,n
         do 50 j=1,m
            read (31,*) nv(1,i,j)
            read (32,*) nv(2,i,j)
            read (33,*) nv(3,i,j)
 50      continue
 100  continue
      do 110 j=3,m-2
         do 105 i=1,n-1
            d=dsqrt(2d0*(1d0+nv(1,i,j)*nv(1,i+1,j)+
     $           nv(2,i,j)*nv(2,i+1,j)+
     $           nv(3,i,j)*nv(3,i+1,j)))
            nv1(1,2*i-1,j)=nv(1,i,j)
            nv1(1,2*i,j)=(nv(1,i,j)+nv(1,i+1,j))/d         
            nv1(2,2*i-1,j)=nv(2,i,j)
            nv1(2,2*i,j)=(nv(2,i,j)+nv(2,i+1,j))/d         
            nv1(3,2*i-1,j)=nv(3,i,j)
            nv1(3,2*i,j)=(nv(3,i,j)+nv(3,i+1,j))/d         
 105     continue
         nv1(1,2*n-1,j)=nv(1,n,j)
         nv1(2,2*n-1,j)=nv(2,n,j)
         nv1(3,2*n-1,j)=nv(3,n,j)
 110  continue
      do 120 i=1,ilimit-1
         d=dsqrt(2d0*(1d0+nv(1,i,2)*nv(1,i+1,2)+
     $        nv(2,i,2)*nv(2,i+1,2)+
     $        nv(3,i,2)*nv(3,i+1,2)))
         nv1(1,2*i-1,2)=nv(1,i,2)
         nv1(1,2*i,2)=(nv(1,i,2)+nv(1,i+1,2))/d         
         nv1(2,2*i-1,2)=nv(2,i,2)
         nv1(2,2*i,2)=(nv(2,i,2)+nv(2,i+1,2))/d         
         nv1(3,2*i-1,2)=nv(3,i,2)
         nv1(3,2*i,2)=(nv(3,i,2)+nv(3,i+1,2))/d         
         d=dsqrt(2d0*(1d0+nv(1,i,m-1)*nv(1,i+1,m-1)+
     $        nv(2,i,m-1)*nv(2,i+1,m-1)+
     $        nv(3,i,m-1)*nv(3,i+1,m-1)))
         nv1(1,2*i-1,m-1)=nv(1,i,m-1)
         nv1(1,2*i,m-1)=(nv(1,i,m-1)+nv(1,i+1,m-1))/d         
         nv1(2,2*i-1,m-1)=nv(2,i,m-1)
         nv1(2,2*i,m-1)=(nv(2,i,m-1)+nv(2,i+1,m-1))/d         
         nv1(3,2*i-1,m-1)=nv(3,i,m-1)
         nv1(3,2*i,m-1)=(nv(3,i,m-1)+nv(3,i+1,m-1))/d         
 120  continue
      nv1(1,2*ilimit-1,2)=nv(1,ilimit,2)   
      nv1(2,2*ilimit-1,2)=nv(2,ilimit,2)   
      nv1(3,2*ilimit-1,2)=nv(3,ilimit,2)   
      nv1(1,2*ilimit-1,m-1)=nv(1,ilimit,m-1)   
      nv1(2,2*ilimit-1,m-1)=nv(2,ilimit,m-1)   
      nv1(3,2*ilimit-1,m-1)=nv(3,ilimit,m-1)   
      do 130 i=2*ilimit,nuusi
         nv1(1,i,2)=0d0
         nv1(2,i,2)=1d0
         nv1(3,i,2)=0d0
         nv1(1,i,m-1)=0d0
         nv1(2,i,m-1)=-1d0
         nv1(3,i,m-1)=0d0
 130  continue
      do 200 i=1,nuusi
         do 150 j=2,m-2
            d=dsqrt(2d0*(1d0+nv1(1,i,j)*nv1(1,i,j+1)+
     $           nv1(2,i,j)*nv1(2,i,j+1)+
     $           nv1(3,i,j)*nv1(3,i,j+1)))
            nv2(1,i,2*(j-1))=nv1(1,i,j)
            nv2(1,i,1+2*(j-1))=(nv1(1,i,j)+nv1(1,i,j+1))/d         
            nv2(2,i,2*(j-1))=nv1(2,i,j)
            nv2(2,i,1+2*(j-1))=(nv1(2,i,j)+nv1(2,i,j+1))/d         
            nv2(3,i,2*(j-1))=nv1(3,i,j)
            nv2(3,i,1+2*(j-1))=(nv1(3,i,j)+nv1(3,i,j+1))/d
 150     continue
         nv2(1,i,2*(m-2))=nv1(1,i,m-1)
         nv2(2,i,2*(m-2))=nv1(2,i,m-1)
         nv2(3,i,2*(m-2))=nv1(3,i,m-1)
         nv2(1,i,1)=nv2(1,i,2*(m-2))
         nv2(2,i,1)=nv2(2,i,2*(m-2))
         nv2(3,i,1)=nv2(3,i,2*(m-2))
         nv2(1,i,2*m-3)=nv2(1,i,2)
         nv2(2,i,2*m-3)=nv2(2,i,2)
         nv2(3,i,2*m-3)=nv2(3,i,2)
 200  continue
      rewind 31
      rewind 32
      rewind 33
      do 400 i=1,nuusi
         do 300 j=1,2*m-3
            write (31,*) nv2(1,i,j)         
            write (32,*) nv2(2,i,j)         
            write (33,*) nv2(3,i,j)
 300     continue
 400  continue
      end
