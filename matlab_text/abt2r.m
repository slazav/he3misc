function r = abt2r(a,b,t)
% convert alpha, beta, theta angles (deg) to the rotation matrix

  n=[sin(b)*cos(a) sin(b)*sin(a) n=cos(b)];

  nn = [ n(1)^2    n(2)*n(1) n(3)*n(1)
         n(1)*n(2) n(2)^2    n(3)*n(2)
         n(1)*n(3) n(2)*n(3) n(3)^2 ];

  en = [    0   n(3) -n(2)
         -n(3)    0   n(1)
          n(2) -n(1)    0];

  dd = [1 0 0; 0 1 0; 0 0 1];

  r = dd*cos(t) + nn*(1-cos(t)) - en*sin(t);
  % or +en ?
end

