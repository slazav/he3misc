function [e1 e2] = en_gr1(a,b, dax,day, dbx,dby, dx,dy)

  n(1) = sin(b)*cos(a);
  n(2) = sin(b)*sin(a);
  n(3) = cos(b);

  ga=[dax/dx day/dy 0];
  gb=[dbx/dx dby/dy 0];

  gn(1,:) = cos(b)*cos(a)*gb - sin(b)*sin(a)*ga;
  gn(2,:) = cos(b)*sin(a)*gb + sin(b)*cos(a)*ga;
  gn(3,:) = -sin(b)*gb;

  e1 = 0;
  e2 = 0;
  u=0;

  c1 = 25/16; % (1-ct)^2
  c2 = 5/8*sqrt(15); % 2*st*(1-ct)
  c3 = 15/16; % st^2

  u = c1 * ( (n(1)*gn(1,1) + n(2)*gn(1,2))^2 + ...
             (n(1)*gn(2,1) + n(2)*gn(2,2))^2 + ...
             (n(1)*gn(3,1) + n(2)*gn(3,2))^2 ) ...
    + c2 * (n(2)*gn(3,1)*gn(1,1) - n(1)*gn(3,2)*gn(2,2) + ...
            n(3) * (gn(1,2)-gn(2,1))*(gn(2,2)+gn(1,1)) );

  e1 = e1 + c2 * (n(2)*gn(3,1)*gn(2,2) - n(1)*gn(3,2)*gn(1,1));
  e2 = e2 + c2 * (n(2)*gn(3,2)*gn(2,1) - n(1)*gn(3,1)*gn(1,2));

  u = u ...
    + 15/16 * (gn(1,2)^2 + gn(2,1)^2 + gn(3,1)^2 + gn(3,2)^2) ...
    + 25/16 * (gn(1,1)^2 + gn(2,2)^2);

  e1 = e1 ...
    + 25/8 * gn(1,1)*gn(2,2) ...
    - 15/8 * gn(1,2)*gn(2,1);

  e2 = e2 ...
    - 15/8 * gn(1,1)*gn(2,2) ...
    + 25/8 * gn(1,2)*gn(2,1);

  e1=e1+u; e2=e2+u;
  return

end

