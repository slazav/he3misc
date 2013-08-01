addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D/')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D/mex/')

addpath('/rota/Analysis/NMRcalc/spinwaves/Lambda_omega/')
addpath('/rota/Analysis/NMRcalc/spinwaves/Lambda_omega/mex/')


Nr=30;Nf=50;
Ntot=Nr*Nf;


Omega=0;Omegav=0;
lo=5;
T=0.15;
f_larmor=826.67*10^3;


nubpar=-1;
Apsi=zeros(Ntot+1,1);

tic;
[textur0,~]=texture_05bar(Nr,[T 0.5 f_larmor*10^-3 0.3 Omega Omegav 5 -1],500,[0.7 10],1);
time=toc

%initial_textur=textur0;
%initial_textur(2:Ntot+1,1:3)=repmat(textur0(2:Nr+1,1:3),Nf,1);

%initial_textur=textur0;
%initial_textur(2:Ntot+1,1:3)=repmat(textur0(2:Nr+1,1:3),Nf,1);
%initial_textur=abs(initial_textur.*(1+(-0.5+rand(size(initial_textur)))*0.2));

%initial_textur(1,2)=rand(1,1)*1000;

initial_textur=1;

tic;
[textur,~]=texture2D(Ntot,[T 0.5 f_larmor*10^-3 0.3 Omega Omegav 5 -1 0 nubpar],500,[0.7 10],initial_textur,Apsi,Nr);

time=toc

figure;hold on
plot(textur(1:Nr+1,1),textur(1:Nr+1,3))
plot(textur(1:Nr+1,1),textur(1:Nr+1,2),'r')
plot(textur0(:,1),textur0(:,3),'--')
plot(textur0(:,1),textur0(:,2),'r--')
xlabel('radius,cm')
ylabel('degrees')

figure;hold on
plot(textur(10:Nr:Ntot+1,4),textur(10:Nr:Ntot+1,3))
plot(textur(10:Nr:Ntot+1,4),textur(10:Nr:Ntot+1,2),'r')
xlabel('phi, rad')
ylabel('degrees')

beta=reshape(textur(2:Ntot+1,3),Nr,Nf);
betaq=beta;
betaq(:,Nf+1)=beta(:,1);
betaq=[textur(1,3);betaq(:)];
beta(2:Nr+1,:)=beta(:,:);
beta(1,:)=textur(1,3);
beta(:,Nf+1)=beta(:,1);


alpha=reshape(textur(2:Ntot+1,2),Nr,Nf);
alphaq=alpha;
alphaq(:,Nf+1)=alpha(:,1);
alphaq=[textur(1,2);alphaq(:)];
alpha(2:Nr+1,:)=alpha(:,:);
alpha(1,:)=textur(1,2);
alpha(:,Nf+1)=alpha(:,1);

r=reshape(textur(2:Ntot+1,1),Nr,Nf);
rq=r;
rq(:,Nf+1)=r(:,1);
rq=[textur(1,1);rq(:)];
r(2:Nr+1,:)=r(:,:);
r(1,:)=textur(1,1);
r(:,Nf+1)=r(:,1);

f=reshape(textur(2:Ntot+1,4),Nr,Nf);
fq=f;
fq(:,Nf+1)=f(:,1);
fq=[textur(1,4);fq(:)];
f(2:Nr+1,:)=f(:,:);
f(1,:)=textur(1,4);
f(:,Nf+1)=f(:,1);

[x,y]=pol2cart(f,r);
[xq,yq]=pol2cart(fq,rq);

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
ncartx=-cosd(-fq*360/(2*pi)).*sind(betaq).*cosd(alphaq) + sind(-fq*360/(2*pi)).*sind(betaq).*sind(alphaq);
ncarty=sind(-fq*360/(2*pi)).*sind(betaq).*cosd(alphaq) + cosd(-fq*360/(2*pi)).*sind(betaq).*sind(alphaq);
quiver(xq,yq,ncartx,ncarty)
pbaspect([1,1,1])
plot(min(xq(:)):(max(xq(:))-min(xq(:)))/100:max(xq(:)),sqrt(max(xq(:))^2-(min(xq(:)):(max(xq(:))-min(xq(:)))/100:max(xq(:))).^2))
plot(min(xq(:)):(max(xq(:))-min(xq(:)))/100:max(xq(:)),-sqrt(max(xq(:))^2-(min(xq(:)):(max(xq(:))-min(xq(:)))/100:max(xq(:))).^2))
