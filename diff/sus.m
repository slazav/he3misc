function s = sus(P,TTc)
  %chi, following D.Einzel JLTP 84 and
  f0a= he3_f0a(P);
  s0 = sus0(P,TTc);
  %calc sus
  s = s0 ./ (1 + f0a * s0);
end

