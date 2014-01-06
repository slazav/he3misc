function e = en_g1(a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz)

  n(1) = sind(b)*cosd(a);
  n(2) = sind(b)*sind(a);
  n(3) = cosd(b);

  ga=[dax/dx day/dy daz/dz];
  gb=[dbx/dx dby/dy dbz/dz];
  gt=[dtx/dx dty/dy dtz/dz];

  gn(1,:) = cosd(b)*cosd(a) * gb - sind(b)*sind(a) * ga;
  gn(2,:) = cosd(b)*sind(a) * gb + sind(b)*cosd(a) * ga;
  gn(3,:) = -sind(b) * gb;

  ee = zeros(3,3,3);
  ee(1,2,3) = 1;
  ee(2,3,1) = 1;
  ee(3,1,2) = 1;
  ee(3,2,1) = -1;
  ee(2,1,3) = -1;
  ee(1,3,2) = -1;

  dd=[1 0 0; 0 1 0; 0 0 1];

  ct=cosd(t);
  st=sind(t);

  dr=zeros(3,3,3);
  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    dr(i,j,k) = dr(i,j,k) + ...
      ((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l))*gn(k,l)  +...
       (st*(n(i)*n(j) - dd(i,j))*dd(i,l) - ct*ee(i,j,l)*n(l))*gt(k);
  end; end; end; end

  e=0;
  for k=1:3; for j=1:3; for i=1:3;
    e = e + dr(i,j,j)*dr(i,k,k);
  end; end; end;

  return

  for k=1:3; for j=1:3; for i=1:3;
    e = e + (1-ct)^2 * n(i)*n(j)*gn(k,i)*gn(k,j);
  end; end; end

  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    e = e - 2*st*(1-ct) * ee(k,j,l)*n(i)*gn(k,i)*gn(l,j);
    e = e - 2*st*(1-ct) * ee(k,j,l)*n(k)*gn(l,j)*gn(i,i);
  end; end; end; end

  for k=1:3; for j=1:3;
    e = e + st^2 * gn(i,j)*gn(i,j);
    e = e + (1-ct)^2 * gn(i,i)*gn(j,j);
    e = e - st^2 * gn(i,j)*gn(j,i);
    e = e - 2*st * n(i)*gn(j,i)*gt(j);
  end; end

  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    e = e - 2*ct*(1-ct) * ee(l,j,k)*n(k)*n(i)*gn(l,i)*gt(j);
    e = e - 2*st^2 * n(k)*n(j)*ee(k,i,l)*gn(l,i)*gt(j);
  end; end; end; end

  for k=1:3; for j=1:3; for i=1:3;
    e = e + 2*st^2 * ee(j,i,k)*gn(k,i)*gt(j);
  end; end; end

  for j=1:3; for i=1:3;
    e = e + (dd(i,j)-n(i)*n(j)) * gt(i)*gt(j);
  end; end




end

