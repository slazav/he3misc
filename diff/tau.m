function [ tau_aver,tauN ] = tau(P,TTc)
%thermal average of quasiparticle lifetime at given energy at low temperature limit
%following D.Einzel JLTP 84, and D. Einzel JLTP 54
kB=8.6173*10^(-5); %eV/K

%parameters
tauN_Tc_vsP=[0.5 0.12 0.054 0.04]*10^(-6);% D.Einzel JLTP 84
PtauN=[ 0 10 20 30];

tauN0=interp1(PtauN,tauN_Tc_vsP,P,'spline');

[T Tc]=TTc_to_T(P,TTc);
addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
sgap=trivial(TTc,P)*kB*Tc;



%calculation
tauN=tauN0*(Tc/T)^2;

phimax=(cosh(sgap/(2*kB*T)))^(-2);
steps=2;
phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);

while phimax<10^4*phi_steps_of_gap
    steps=steps*1.2;
 phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);
end

%calc integral
tau_aver=tauN./(2/(4*kB*T)*quad(@(ksi) integrand(ksi,sgap,kB,T,P),0,sgap*steps,10^-8));%yoshida0 omitted in Itau
%2*quad(@(ksi) integrand(ksi,sgap,kB,T,P),0,sgap*steps,10^-8)

%ksi_test=0:sgap*steps/100:sgap*steps;figure;plot(ksi_test,integrand(ksi_test,sgap,kB,T,P));
end

function integr=integrand(ksi,sgap,kB,T,P)
    Ek=sqrt(ksi.^2+sgap^2);
    phi=(cosh(Ek./(2*kB*T))).^(-2);
    integr=phi.*Itau(ksi,P,T,sgap);
end


function Itau=Itau(ksi,P,T,sgap)
   %implementation of energy dependence of bogoliubov quasiparticles at low temperatures T<<Tc,
   %following D.Einzel JLTP 84, and D. Einzel JLTP 54
    
   %parameters
    gamma0=0.1; %only weak press. dep.
    delta0=0.3; %only weak press. dep.
    w0=1-2*gamma0/3+delta0;
        
    kB=8.6173*10^(-5); %eV/K
    
      
    Izero=3*sgap/(2*pi*kB*T);  %omit yoshida0, it is cancelled in the average life time calculation
    x=ksi./sqrt(2*kB*T*sgap);
    Itau=Izero*(w0+kB*T/sgap*(3/4*(1+x.^2)*w0-(1+2*x.^2)*(gamma0/3+delta0)));   %%TESTING!!!!!!!!!!!!!!!!!!!!!!!!
end

