format long;
N=250; %number of steps: NOT TO BE VARIED IN TEXTURE CALCULATION
R=0.3*10^(-2);
h=R/(N+1);  % h=step size
r= (0:h:((N+1)*h))';

chi=1.437*10^(-8);

T=[0.14 0.18 0.22 0.26]; %    T / Tc
lambda_omega=5.5;

%To frequences
for k=1:length(T)
f_larmor=826.67*10^3; %Herz
f_L=sqrt(14.46/16.8075*(1-T(k).^2).*(44.2121*T(k).^6-64.5411*T(k).^4+16.9909*T(k).^2+16.862)*1000)*10^3; %Herz
C(k)=f_L^2/(2*f_larmor);
end
%Caluculation paratemeters

omega=[0.5];
omegav=[0.5];

steps_temperature=length(T);
kmax=25;

steps=100;
M_max=60*10^(-3);

M=0:M_max/steps:(steps-1)/steps*M_max;



FR=zeros(steps,steps_temperature);
ratios=zeros(steps,steps_temperature);
npsi_previous=zeros(N+2,steps,steps_temperature);
textures=zeros(N+2,3,steps,steps_temperature);

for s=1:steps_temperature   
    s
    FRtemp=zeros(steps,1); 
    ratios_temp=zeros(steps,1);
    npsi_previous_temp=zeros(N+2,steps);
    previous_textures=zeros(N+2,3,steps);
for n=1:steps
    n
    averageFR=0;
    k=1;
        
    %Initial value for npsi
    
    if n==1
        [E,Eigenvectors]=energies_M05bar(N,T(s),f_larmor/1000,lambda_omega,omega,omegav,chi,zeros(N+2,1),1);
        psi=Eigenvectors(:,1);  
        psi(2:N+1)=psi;
        psi(1)=psi(2);
        psi(N+2)=0;
        [textur,~]=textureM05bar(N+1,[T(s) 0.5 f_larmor/1000 0.3 omega omegav 5.5 0.9 chi],500,[0.7 10],1,zeros(N+2,1));
        initial_texture=textur;

    else       
       initial_texture=previous_textures(:,:,n-1);      
       psi=npsi_previous_temp(:,n-1);
    end
    
    ratio=1;    
    previous_vectors=zeros(N+2,kmax+2);
    previousFR=zeros(kmax+2,1);
    previous_ratio=1;
    
    while  k<=kmax
        

        %Normalization over volume WITHOUT 2*pi*h
        norm=trapz(r,r.*(psi).^2);
        npsi=psi/sqrt(norm);
      
        %using averaged npsi from all rounds
        previous_vectors(:,k)=abs(npsi); 
        
        if k>1
            averaged=sum(transpose(previous_vectors(:,1:k)))./k;
            averaged=averaged';
            [A res]=findA(N,M(n),averaged);
            [E,Eigenvectors]=energies_M05bar(N,T(s),f_larmor/1000,lambda_omega,omega,omegav,chi,A*averaged,initial_texture);
        else
            [A res]=findA(N,M(n),npsi);
            [E,Eigenvectors]=energies_M05bar(N,T(s),f_larmor/1000,lambda_omega,omega,omegav,chi,A*npsi,initial_texture);
        end
        psi=Eigenvectors(:,1);  
        psi(2:N+1)=psi;
        psi(1)=psi(2);
        psi(N+2)=0;        
        
        FR1=E(1)*C(s); 
        previousFR(k)=FR1; 
        
        if k>1
            averageFR=sum(previousFR(1:k))/k;
            ratio=abs((FR1-averageFR)/FR1);
        end
                
        previous_ratio=ratio;        
        k=k+1
    end
    
    FRtemp(n,1)=averageFR;     
    ratios_temp(n,1)=ratio;
    npsi_previous_temp(:,n)=sum(transpose(previous_vectors(:,1:k-1)))./(k-1);
    [textur1,~]=textureM05bar(N+1,[T(s) 0.5 f_larmor/1000 0.3 omega omegav 5.5 0.9 chi],500,[0.7 10],initial_texture,A*averaged);
    previous_textures(:,:,n)=textur1;

end
FR(:,s)=FRtemp;
ratios(:,s)=ratios_temp;
npsi_previous(:,:,s)=npsi_previous_temp;
textures(:,:,:,s)=previous_textures;
end



%plot results
figure
hold on
for s=1:steps_temperature
    plot(FR(:,s).*10^-3,M)
end





