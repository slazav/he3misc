function [ yos,steps ] = yosida(P,TTc,n)
%IMPLEMENTATION OF YOSIDA FUNCTIONS, FOLLOWING D.Einzel JLTP 84

%parameters
kB=8.6173*10^(-5); %eV/K
[T Tc]=TTc_to_T(P,TTc);

addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
sgap=trivial(TTc,P)*kB*Tc;




%integral
%find relevant energy scale
phimax=(cosh(sgap/(2*kB*T)))^(-2);
steps=2;
phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);

while phimax<10^4*phi_steps_of_gap
    steps=steps*1.2;
 phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);
end

%calc integral
integ=2*quad(@(ksi) integrand(ksi,sgap,kB*T,n),0,sgap*steps,10^-8);
yos=integ/(4*kB*T);


end

function integr = integrand(ksi,sgap,kBT,n)
    Ek=sqrt(ksi.^2+sgap^2);
    phi=(cosh(Ek./(2*kBT))).^(-2);
    integr=abs(ksi./Ek).^n.*phi;
end