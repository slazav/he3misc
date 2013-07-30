%CLEAR (delete) frot.25 frot.27  BEFORE EXECUTING and boot matlab
Nr=30;Nf=50;
Ntot=Nr*Nf;

if(1==0)
load('testing_DONOTDELETE.dat')
testing=testing_DONOTDELETE;
gr_2Dt=testing(length(testing)-1501:length(testing)-1,:);
initial_texture=gr_2Dt(1:31,1:3);
initial_texture(:,2:3)=initial_texture(:,2:3)*180/pi;
initial_texture(20:21,2:3)=initial_texture(30:31,2:3);
initial_texture(10:11,2:3)=initial_texture(1:2,2:3);

else
initial_texture=zeros(Nr+1,3);
initial_texture(Nr+1:-1:1,2)=170:-170/Nr:0;
initial_texture(Nr+1:-1:1,3)=0:64/Nr:64;
end



initial_textur_2D=initial_texture;
initial_textur_2D(2:Ntot+1,1:3)=repmat(initial_texture(2:Nr+1,1:3),Nf,1);

addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/Original_Kopus/kopu_1D/lib')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/Original_Kopus/kopu_1D/mex')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/mex/')



Omega=1.4;Omegav=0;
lo=-1;
T=0.58;
f_larmor=910*10^3;  %corresponds to 28mT
P=29.3;

%nubpar=-1;
nubpar=100*sqrt((1-T^4)*(9.00301-19.927*T^4+15.3442*T^6)); %khz
Apsi=zeros(Ntot+1,1);

% tic;
% [textur0,~]=texture_05bar(Nr,[T P f_larmor*10^-3 0.3 Omega Omegav 5 -1],500,[0.7 10],1);
% time=toc
% 

[textur,~]=texture_testing(Nr,[T P f_larmor*10^-3 0.3 Omega Omegav lo -1],500,[0.7 10],initial_texture);
%%
[textur_2D,~]=texture2D(Ntot,[T P f_larmor*10^-3 0.3 Omega Omegav lo -1 0 nubpar],500,[0.7 10],initial_textur_2D,Apsi,Nr);

%% 
%files fort.25 and fort.27 NEED TO BE HAND MODIFIED AT THIS POINT


gr_1D=load('fort.25.dat');
gr_2D=load('fort.27.dat');

[gr_1D(:,[1,4,5]),gr_2D(1:Nr+1,4:5)*Nf/(2*pi)]
figure;
subplot(2,1,2)
hold on
plot(gr_2D(1:Nr+1,4)*Nf/(2*pi),'r')
plot(gr_2D(1:Nr+1,5)*Nf/(2*pi))
plot(gr_1D(:,[4]),'r--')
plot(gr_1D(:,[5]),'--')

subplot(2,1,1);hold on
plot(gr_2D(1:Nr+1,2),'r')
plot(gr_2D(1:Nr+1,3))
plot(gr_1D(:,[2]),'r--')
plot(gr_1D(:,[3]),'--')


%% 
%this is to create plot v1


gr_1D=load('fort.25_v1.dat');
gr_2D=load('fort.27_v1.dat');

[gr_1D(:,[1,4,5]),gr_2D(1:Nr+1,4:5)*Nf/(2*pi)]
figure;
subplot(2,1,2)
hold on
plot(gr_2D(1:Nr+1,4)*Nf/(2*pi),'r')
plot(gr_2D(1:Nr+1,5)*Nf/(2*pi))
plot(gr_1D(:,[4]),'r--')
plot(gr_1D(:,[5]),'--')

subplot(2,1,1);hold on
plot(gr_2D(1:Nr+1,2),'r')
plot(gr_2D(1:Nr+1,3))
plot(gr_1D(:,[2]),'r--')
plot(gr_1D(:,[3]),'--')


%% 
%this is to create plot v2


gr_1D=load('fort.25_v2.dat');
gr_2D=load('fort.27_v2.dat');

[gr_1D(:,[1,4,5]),gr_2D(1:Nr+1,4:5)*Nf/(2*pi)]
figure;
subplot(2,1,2)
hold on
plot(gr_2D(1:Nr+1,4)*Nf/(2*pi),'r')
plot(gr_2D(1:Nr+1,5)*Nf/(2*pi))
plot(gr_1D(:,[4]),'r--')
plot(gr_1D(:,[5]),'--')

subplot(2,1,1);hold on
plot(gr_2D(1:Nr+1,2),'r')
plot(gr_2D(1:Nr+1,3))
plot(gr_1D(:,[2]),'r--')
plot(gr_1D(:,[3]),'--')