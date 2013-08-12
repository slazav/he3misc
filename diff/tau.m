function [ tau_aver,tauN ] = tau(P,TTc)
  %thermal average of quasiparticle lifetime at given energy at low temperature limit
  %following D.Einzel JLTP 84, and D. Einzel JLTP 54

  %parameters
  tauN_Tc_vsP=[0.5 0.12 0.054 0.04]*10^(-6);% D.Einzel JLTP 84
  PtauN=[ 0 10 20 30];
  tauN0=interp1(PtauN,tauN_Tc_vsP,P,'spline');

  gap=he3_trivgap(TTc,P);

  %calculation
  tauN=tauN0/TTc^2;

  %calc integral
  tau_aver=tauN./(2/(4*TTc)*quad(@(x) integrand(x,gap,TTc,P),0,1,10^-8));
  %yoshida0 omitted in Itau
end

function integr=integrand(x,gap,TTc,P)

   % Itau is implementation of energy dependence of bogoliubov
   % quasiparticles at low temperatures T<<Tc,
   % following D.Einzel JLTP 84, and D. Einzel JLTP 54

   %parameters
    gamma0=0.1; %only weak press. dep.
    delta0=0.3; %only weak press. dep.
    w0 = 1 - 2/3*gamma0 + delta0;

    ksi = atanh(x) *(2*TTc);

    Izero=3*gap/(2*pi*TTc);  %omit yoshida0, it is cancelled in the average life time calculation
    xx=ksi./sqrt(2*TTc*gap);
    Itau=Izero*(w0+TTc/gap*(3/4*(1+xx.^2)*w0-(1+2*xx.^2)*(gamma0/3+delta0)));   %%TESTING!!!!!!!!!!!!!!!!!!!!!!!!

    Ek=sqrt(ksi.^2+gap^2);
    phi=(cosh(Ek./(2*TTc))).^(-2);

    integr = phi .* Itau * (2*TTc) .* cosh(ksi/(2*TTc)).^2;
end


