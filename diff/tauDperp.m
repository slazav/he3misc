function [ tauDperp ] = tauDperp( TTc,P )
  [ tau_aver,~] = tau(P,TTc);
  %following D.Einzel JLTP 84, and D. Einzel JLTP 54
  lambda=he3_scatt_l1a(P);
  gap = he3_trivgap(TTc,P);
  y0=he3_yosida(TTc, gap, 0);
  y2=he3_yosida(TTc, gap, 2);
  tauDperp=tau_aver/(1-lambda*(4/5*y0+1/5*y2)/y0);
end
