function e = en_d_0(a,b,t)
  r = abt2r(a, b, t);
  e = 1/5 * (trace(r^2) + trace(r)^2);
end

