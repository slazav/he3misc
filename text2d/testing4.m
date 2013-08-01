%CLEAR (delete) frot.25 BEFORE EXECUTING

load('testing_DONOTDELETE.dat')
testing=testing_DONOTDELETE;
gr_2D=testing(length(testing)-1501:length(testing)-1,:);
initial_texture=gr_2D(1:31,1:3);
initial_texture(:,2:3)=initial_texture(:,2:3)*180/pi;

addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/Original_Kopus/kopu_1D/lib')
addpath('/rota/Analysis/NMRcalc/spinwaves/TextLib2D_tilt/Original_Kopus/kopu_1D/mex')


Nr=30;Nf=50;
Ntot=Nr*Nf;


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

[textur0,spec]=texture_testing(Nr,[T P f_larmor*10^-3 0.3 Omega Omegav lo -1],500,[0.7 10],initial_texture);
%%
%fort.25 needs to be HAND MODIFIED AT THIS STAGE INTO fort.25.dat by
%removing additional lines

gr_1D=load('fort.25_DNR.dat');


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
