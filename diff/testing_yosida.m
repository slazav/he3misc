
clear all


TTc=0.1:0.01:0.8;

P=0.5;

kB=8.6173*10^(-5); %eV/K
[T Tc]=TTc_to_T(P,TTc);

addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
for k=1:length(TTc)
sgap(k)=trivial(TTc(k),P)*kB*Tc;
gap_kBT_vec(k)=gap_kBT(TTc(k),P);
end


n=0;
for k=1:length(TTc)
   yos(k)=yosida(P,TTc(k),n); 
end
plot(TTc,yos)
yos_lowT=2*gamma((n+1)/2)*(2*kB*T./sgap).^((n-1)/2) .*exp(-sgap./(kB*T));
figure;loglog(gap_kBT_vec,yos./yos_lowT)
hold on
loglog(gap_kBT_vec,gap_kBT_vec.^2)

figure;hold on
plot(TTc,yos)
plot(TTc,yos_lowT,'--')



n=1;
for k=1:length(TTc)
   yos(k)=yosida(P,TTc(k),n); 
end
plot(TTc,yos)
yos_lowT=2*gamma((n+1)/2)*(2*kB*T./sgap).^((n-1)/2) .*exp(-sgap./(kB*T));
plot(TTc,yos_lowT,'--')


n=2;
for k=1:length(TTc)
   yos(k)=yosida(P,TTc(k),n); 
end
plot(TTc,yos)
yos_lowT=2*gamma((n+1)/2)*(2*kB*T./sgap).^((n-1)/2) .*exp(-sgap./(kB*T));
plot(TTc,yos_lowT,'--')


n=3;
for k=1:length(TTc)
   yos(k)=yosida(P,TTc(k),n); 
end
plot(TTc,yos)
yos_lowT=2*gamma((n+1)/2)*(2*kB*T./sgap).^((n-1)/2) .*exp(-sgap./(kB*T));
plot(TTc,yos_lowT,'--')


n=4;
for k=1:length(TTc)
   yos(k)=yosida(P,TTc(k),n); 
end
plot(TTc,yos)
yos_lowT=2*gamma((n+1)/2)*(2*kB*T./sgap).^((n-1)/2) .*exp(-sgap./(kB*T));
plot(TTc,yos_lowT,'--')

