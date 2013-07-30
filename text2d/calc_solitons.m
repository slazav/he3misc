addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib_Slava1/')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib_Slava1/mex/')

addpath('/rota/Analysis/NMRcalc/spinwaves/Lambda_omega/')
addpath('/rota/Analysis/NMRcalc/spinwaves/Lambda_omega/mex/')


Nr=50;Nf=50;
Ntot=Nr*Nf;


Omega=2;Omegav=0.5;
lo=5;
T=0.7;
f_larmor0=826.67*10^3;
f_larmor=823*10^3;

nub=100*sqrt((1-T^4)*(9.00301-19.927*T^4+15.3442*T^6)); %khz
%nubpar=-1;
Apsi=zeros(Ntot+1,1);

tic;
[textur0,~]=texture_05bar(Nr,[T 0.5 f_larmor0*10^-3 0.3 Omega Omegav 5 -1],500,[0.7 10],1);
time=toc

initial_textur=textur0;
initial_textur(2:Ntot+1,1:3)=repmat(textur0(2:Nr+1,1:3),Nf,1);
initial_textur(2:Ntot+1,2)=-initial_textur(2:Ntot+1,2);
%initial_textur=initial_textur.*(1+(-0.5+rand(size(initial_textur)))*0.05);
%initial_textur=1;
tic;
[textur,~]=texture2D_05bar(Ntot,[T 29 f_larmor*10^-3 0.29 Omega Omegav 5 -1 0 nub],500,[0.7 10],initial_textur,Apsi,Nr);

time=toc

figure;hold on
plot(textur(1:Nr+1,1),textur(1:Nr+1,3))
plot(textur(1:Nr+1,1),textur(1:Nr+1,2),'r')
plot(textur0(:,1),textur0(:,3),'--')
plot(textur0(:,1),textur0(:,2),'r--')
xlabel('radius,cm')
ylabel('degrees')

beta=reshape(textur(2:Ntot+1,3),Nr,Nf);
beta(2:Nr+1,:)=beta(:,:);
beta(1,:)=textur(1,3);
beta(:,Nf+1)=beta(:,1);
alpha=reshape(textur(2:Ntot+1,2),Nr,Nf);
alpha(2:Nr+1,:)=alpha(:,:);
alpha(1,:)=textur(1,2);
alpha(:,Nf+1)=alpha(:,1);
r=reshape(textur(2:Ntot+1,1),Nr,Nf);
r(2:Nr+1,:)=r(:,:);
r(1,:)=textur(1,1);
r(:,Nf+1)=r(:,1);
f=reshape(textur(2:Ntot+1,4),Nr,Nf);
f(2:Nr+1,:)=f(:,:);
f(1,:)=textur(1,4);
f(:,Nf+1)=f(:,1);
[x,y]=pol2cart(f,r);

figure;subplot(2,1,1)
surf(x,y,beta)
zlabel('beta, degrees')
xlabel('x,cm')
ylabel('y,cm')
subplot(2,1,2);
surf(x,y,alpha)
zlabel('alpha, degrees')
xlabel('x,cm')
ylabel('y,cm')

figure;hold on
quiver(x,y,sind(beta)*cosd(alpha+f*360/(2*pi)),sind(beta)*sind(alpha+f*360/(2*pi)))
pbaspect([1,1,1])
plot(min(x(:)):(max(x(:))-min(x(:)))/100:max(x(:)),sqrt(max(x(:))^2-(min(x(:)):(max(x(:))-min(x(:)))/100:max(x(:))).^2))
plot(min(x(:)):(max(x(:))-min(x(:)))/100:max(x(:)),-sqrt(max(x(:))^2-(min(x(:)):(max(x(:))-min(x(:)))/100:max(x(:))).^2))
