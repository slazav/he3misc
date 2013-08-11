function [ T,Tc] = TTc_to_T(P,TTc )
  Tc=he3_tc(P)/1000;
  T=TTc*Tc;
end

