function [ tau_aver,tauN ] = tau(P,TTc)
  %thermal average of quasiparticle lifetime at given energy at low temperature limit
  %following D.Einzel JLTP 84, and D. Einzel JLTP 54

  %parameters
  for i = 1:length(TTc)
    tauN=he3_tau_n0tc(P) ./ TTc(i).^2;
    gap=he3_trivgap(TTc(i),P);
    %calc integral
    tau_aver(i)=tauN./(2./(4*TTc(i))*quad(@(x) integrand(x,gap,TTc(i),P),0,1,10^-8));
  end
end

function integr=integrand(x,gap,TTc,P)

   % Itau is implementation of energy dependence of bogoliubov
   % quasiparticles at low temperatures T<<Tc,
   % following D.Einzel JLTP 84, and D. Einzel JLTP 54

   %parameters
    gamma0=0.1; %only weak press. dep.
    delta0=0.3; %only weak press. dep.
    w0 = 1 - 2/3*gamma0 + delta0;

    ksi = atanh(x)*(2*TTc);

    Izero=3*gap/(2*pi*TTc);  %omit yoshida0, it is cancelled in the average life time calculation
    xx=ksi./sqrt(2*TTc*gap);
    Itau=Izero*(w0+TTc/gap*(3/4*(1+xx.^2)*w0-(1+2*xx.^2)*(gamma0/3+delta0)));   %%TESTING!!!!!!!!!!!!!!!!!!!!!!!!

    Ek=sqrt(ksi.^2+gap^2);
    phi=(cosh(Ek./(2*TTc))).^(-2);

    integr = phi .* Itau * (2*TTc) .* cosh(ksi/(2*TTc)).^2;
end


