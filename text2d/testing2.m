
Nr=100;
Nf=100;
dx=1/Nr;
df=2*pi/Nf;

sp=(3 + sqrt(3))/6;
sm=(3 - sqrt(3))/6;

      spm(:,1,1) = [ sp*sp,sp*sm,sp*sm,sm*sm ]; %!to first interp. point from 4 corner points
      spm(:,2,1) = [ sp*sm,sp*sp,sm*sm,sp*sm ]; %!second interp. point
      spm(:,1,2) = [ sp*sm,sm*sm,sp*sp,sp*sm ]; %!third  interp. point
      spm(:,2,2) = [ sm*sm,sm*sp,sm*sp,sp*sp ]; %!fourth interp. point1
      
      ind(:,1)= [ 0,1,0,1 ];
      ind(:,2)=[ 0,0,1,1 ];%!f-direction
        


load('testing.dat')

%%

alpha_mat=reshape(testing(1:4:Nr*(Nf)*4,6),Nf,Nr);
alpha_mat=alpha_mat';
alpha_mat(2:Nr+1,:)=alpha_mat;
alpha_mat(1,:)= acos(-1/sqrt(5))/4;
alpha_mat(:,Nf+1)=alpha_mat(:,1);

r=reshape(testing(1:4:Nr*(Nf)*4,1),Nf,Nr);
r=r';
r(2:Nr+1,:)=r;
r(1,:)= 0;
r(:,Nf+1)=r(:,1);
r=r/Nr;

en=0;
en2=0;
for i=1:Nr
     rpm(2)=(i-1+sp)*dx;
     rpm(1)=(i-1+sm)*dx;
for ii=1:Nf        
    apm=0;
    for k=1:2       
    for kk=1:2
    for kkk=1:4
     in1=i+ind(kkk,1);
     in2=ii+ind(kkk,2);
     
     apm = apm + spm(kkk,k,kk)*alpha_mat(in1,in2);
      
    end
    en = en + rpm(k)*(sin(apm))^2*0.25*dx*df;
    
    %e = e + rpm(k)*0.25*dx*df;
    end
    end
    
end 
end
en
