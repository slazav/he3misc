function e = en_d_0(a,b,t)
%  e=0;
%  for i=1:3; for j=1:3;
%    e = e + r(i,i)*r(j,j) + r(i,j)*r(i,j) + r(i,j)*r(j,i);
%  end; end

  r = abt2r(a, b, t);
  e = (trace(r^2) + trace(r)^2) + 3;

end

