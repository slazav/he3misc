format long;
N=200; %number of steps: NOT TO BE VARIED IN TEXTURE CALCULATION
R=0.3*10^(-2);
h=R/(N+1);  % h=step size
r= (0:h:((N+1)*h))';

chi=1.437*10^(-8);

%Caluculation paratemeters

Temp=[0.14]; %    T / Tc

lambda_omega=-1;
lambdaHV=-1;

omegas=[0.5 2 20 100];
%omegas=1;
omegavs=zeros(length(omegas),1);%omegas;%[0];% 0 1];

steps_omega=length(omegas);
kmax=20;

steps=1500;
M_max=1500*10^(-3);

M=0:M_max/steps:(steps-1)/steps*M_max;



FR=zeros(steps,steps_omega);
ratios=zeros(steps,steps_omega);
npsi_previous=zeros(N+2,steps,steps_omega);
Aiter=zeros(steps,steps_omega);
textures=zeros(N+2,3,steps,steps_omega);
poolsize=8;
if matlabpool('size')==0    
    matlabpool('local', num2str(poolsize))
elseif matlabpool('size')<poolsize
    matlabpool close
    matlabpool('local', num2str(poolsize))
end

parfor s=1:steps_omega   
    s;
    T=Temp;
    
    %To frequences
    f_larmor=826.67*10^3; %Herz
    f_L=sqrt(14.46/16.8075*(1-T.^2).*(44.2121*T.^6-64.5411*T.^4+16.9909*T.^2+16.862)*1000)*10^3; %Herz
    C=f_L^2/(2*f_larmor);
    
    [textur_calc,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 0.1 0.1 lambda_omega lambdaHV chi],500,[0.7 10],1,zeros(N+2,1));
    omega=omegas(s);
    omegav=omegavs(s);
    FRtemp=zeros(steps,1); 
    ratios_temp=zeros(steps,1);
    npsi_previous_temp=zeros(N+2,steps);
    previous_textures=zeros(N+2,3,steps);
    Aiter_temp=zeros(steps,1);
for n=1:steps
    n
    averageFR=0;
    k=1;
        
    %Initial value for npsi
    
    if n==1
        [E,Eigenvectors]=energies_M05bar_largeangles(N,T,f_larmor,omega,omegav,lambdaHV,chi,zeros(N+2,1),textur_calc);
        psi=Eigenvectors(:,1);  
        psi(2:N+1)=psi;
        psi(1)=psi(2);
        psi(N+2)=0;
        [textur,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 omega omegav lambda_omega lambdaHV chi],500,[0.7 10],textur_calc,zeros(N+2,1));
        for dummy=1:5
        [textur,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 omega omegav lambda_omega lambdaHV chi],500,[0.7 10],textur,zeros(N+2,1));
        end
        initial_texture=textur;

    elseif n>4
        initial_texture=previous_textures(:,:,n-1);      
        psi=npsi_previous_temp(:,n-1);
        for tt=1:3
            initial_texture=initial_texture+previous_textures(:,:,n-1-tt);      
            psi=psi+npsi_previous_temp(:,n-1-tt);
        end
        initial_texture=initial_texture./4;      
        psi=psi./4;
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
        previous_vectors(:,k)=npsi./npsi(1); 
        
        if k>1
            averaged=sum(transpose(previous_vectors(:,1:k)))./k;
            averaged=averaged';
            [A res]=findA(N,M(n),averaged);
            [E,Eigenvectors]=energies_M05bar_largeangles(N,T,f_larmor,omega,omegav,lambdaHV,chi,A*averaged,initial_texture);
        else
            [A res]=findA(N,M(n),npsi);
            [E,Eigenvectors]=energies_M05bar_largeangles(N,T,f_larmor,omega,omegav,lambdaHV,chi,A*npsi,initial_texture);
        end
        psi=Eigenvectors(:,1);  
        psi(2:N+1)=psi;
        psi(1)=psi(2);
        psi(N+2)=psi(N+1);        
        
        FR1=E(1)*C; 
        previousFR(k)=FR1; 
        
        if k>1
            averageFR=sum(previousFR(1:k))/k;
            ratio=abs((FR1-averageFR)/FR1);
        end
                
        previous_ratio=ratio        
        k=k+1
    end
    
    FRtemp(n,1)=averageFR;
    Aiter_temp(n)=A;
    ratios_temp(n,1)=ratio;
    npsi_previous_temp(:,n)=sum(transpose(previous_vectors(:,1:k-1)))./(k-1);
    [textur1,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 omega omegav lambda_omega lambdaHV chi],500,[0.7 10],initial_texture,A*averaged);
    for dummy=1:5
        [textur1,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 omega omegav lambda_omega lambdaHV chi],500,[0.7 10],textur1,A*averaged);
    end
    previous_textures(:,:,n)=textur1;

end
FR(:,s)=FRtemp;
ratios(:,s)=ratios_temp;
npsi_previous(:,:,s)=npsi_previous_temp;
textures(:,:,:,s)=previous_textures;
Aiter(:,s)=Aiter_temp;
end
matlabpool close


%plot results
figure
hold on
for s=1:steps_omega
    plot(FR(:,s).*10^-3,M)
end
xlabel('freq. kHz')
ylabel('M')




