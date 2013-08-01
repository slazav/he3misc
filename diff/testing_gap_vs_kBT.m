function gap_kBT=gap_kBT(TTc,P)

kB=8.6173*10^(-5); %eV/K
[T Tc]=TTc_to_T(P,TTc);

addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
sgap=trivial(TTc,P)*kB*Tc;

gap_kBT=sgap/(kB*T)
end

