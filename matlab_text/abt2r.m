function r = abt2r(a,b,t)
% convert alpha, beta, theta angles (deg) to the rotation matrix

  n=[sind(b)*cosd(a) sind(b)*sind(a) n=cosd(b)];

  nn = [ n(1)^2    n(2)*n(1) n(3)*n(1)
         n(1)*n(2) n(2)^2    n(3)*n(2)
         n(1)*n(3) n(2)*n(3) n(3)^2 ];

  en = [     0 -n(3)  n(2)
          n(3)    0  -n(1)
         -n(2)  n(1)    0];

  d = [1 0 0; 0 1 0; 0 0 1];

  r = d*cosd(t) + nn*(1-cosd(t)) - en * sind(t);
  % or +en ?
end

