Nksi=30;
poolsize=8;
if matlabpool('size')==0    
    matlabpool('local', num2str(poolsize))
elseif matlabpool('size')<poolsize
    matlabpool close
    matlabpool('local', num2str(poolsize))
end

TTc=0.07:0.01:0.7;

w0=826*10^3*2*pi;
Press=0.5;


kB=8.6173*10^(-5); %eV/K

steps_vec=zeros(size(TTc));
for k=1:length(TTc)   
[T Tc]=TTc_to_T(Press,TTc(k));

addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
sgap(k)=trivial(TTc(k),Press)*kB*Tc;




phimax=(cosh(sgap(k)/(2*kB*T)))^(-2);
steps=2;
phi_steps_of_gap=(cosh(sgap(k)*steps/(2*kB*T)))^(-2);

while phimax<2*10^6*phi_steps_of_gap
    steps=steps*1.2;
 phi_steps_of_gap=(cosh(sgap(k)*steps/(2*kB*T)))^(-2);
end
steps_vec(k)=steps;
end


Diffc=cell(length(TTc),1);
ksi_vec_cell=cell(length(TTc),1);
for k=1:length(TTc)   
    ksi_vec=0:steps_vec(k)*sgap(k)/Nksi:steps_vec(k)*sgap(k);
    ksi_vec_cell{k}=ksi_vec;
    Diffc_help=zeros(length(ksi_vec),1);
parfor kk=1:length(ksi_vec) 
    
    ksi=ksi_vec(kk);
   Diffc_help(kk)=angle_integr_lowtemp(Press,TTc(k),w0,ksi); 
end
Diffc{k}=Diffc_help;
end



%%
figure;hold on

for k=1:length(TTc)
    plot(ksi_vec_cell{k},Diffc{k})
    Ek=sqrt(ksi_vec_cell{k}.^2+sgap(k)^2);
    ksi_E=ksi_vec_cell{k}./Ek;
    ksi_E=ksi_E.^2*10^7;
    plot(ksi_vec_cell{k},ksi_E,'--')
end