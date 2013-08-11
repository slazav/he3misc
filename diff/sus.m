function [ sus] = sus(P,TTc)
  %chi, following D.Einzel JLTP 84 and
  f0a= he3_f0a(P);
  %calc sus
  sus = sus0(P,TTc)/(1 + f0a * sus0(P,TTc));
end

