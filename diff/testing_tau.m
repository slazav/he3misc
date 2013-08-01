clear all
TTc=0.05:0.01:0.99;
P=0.5;
for k=1:length(TTc)
    [tau_vec(k) tauN(k)]=tau(P,TTc(k));
    gap_kBT_vec(k)=gap_kBT(TTc(k),P);

end


figure;hold on
plot(gap_kBT_vec,tauN./tau_vec)


%low temp limit according to D.Einzel JLTP 84

%parameters
    gamma0=0.1; %only weak press. dep.
    delta0=0.3; %only weak press. dep.
    w0=1-2*gamma0/3+delta0;
    kB=8.6173*10^(-5); %eV/K
    [T Tc]=TTc_to_T(P,TTc);
    addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
    
    for k=1:length(TTc)
    sgap(k)=trivial(TTc(k),P)*kB*Tc;
    Izero(k)=3*sgap(k)*yosida(P,TTc(k),0)/(2*pi*kB*T(k));
    %Izero(k)=yosida(P,TTc(k),0)/w0
    end
    
    
    
    tau_lowT=tauN./(Izero*w0);

plot(gap_kBT_vec,tauN./tau_lowT,'--')
%%
figure;loglog(gap_kBT_vec,tau_vec./tau_lowT)
hold on
loglog(gap_kBT_vec,gap_kBT_vec.^2)