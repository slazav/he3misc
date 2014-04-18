function [e1 e2] = en_gr(a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz)

  n(1) = sin(b)*cos(a);
  n(2) = sin(b)*sin(a);
  n(3) = cos(b);

  ga=[dax/dx day/dy daz/dz];
  gb=[dbx/dx dby/dy dbz/dz];
  gt=[dtx/dx dty/dy dtz/dz];

  gn(1,:) = cos(b)*cos(a)*gb - sin(b)*sin(a)*ga;
  gn(2,:) = cos(b)*sin(a)*gb + sin(b)*cos(a)*ga;
  gn(3,:) = -sin(b)*gb;

  ee = zeros(3,3,3);
  ee(1,2,3) = 1;
  ee(2,3,1) = 1;
  ee(3,1,2) = 1;
  ee(3,2,1) = -1;
  ee(2,1,3) = -1;
  ee(1,3,2) = -1;

  dd=[1 0 0; 0 1 0; 0 0 1];

  ct=cos(t);
  st=sin(t);

  %%% variant 1 - best for computations?
%  gr = zeros(3,3,3);
%  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
%    gr(i,j,k) = gr(i,j,k) + ...
%      ((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l)) * gn(l,k) + ...
%      (st*(n(i)*n(j) - dd(i,j))*dd(j,l) - ct*ee(i,j,l)*n(l))*gt(k);
%  end; end; end; end;

%  e1 = 0;
%  e2 = 0;
%  for k=1:3; for j=1:3; for i=1:3;
%    e1 = e1 + gr(i,j,j) * gr(i,k,k);
%    e2 = e2 + gr(i,j,k) * gr(i,k,j);
%  end; end; end;
%  return


%  %%% variant 2
%  e1 = 0;
%  e2 = 0;
%  for l=1:3; for m=1:3; for k=1:3; for j=1:3; for i=1:3;
%    e1 = e1 +...
%      (((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l)) * gn(l,j) + ...
%      (st*(n(i)*n(j) - dd(i,j))*dd(j,l) - ct*ee(i,j,l)*n(l))*gt(j)) *...
%       (((1-ct)*(dd(i,m)*n(k) + dd(k,m)*n(i)) - st*ee(i,k,m)) * gn(m,k) + ...
%      (st*(n(i)*n(k) - dd(i,k))*dd(k,m) - ct*ee(i,k,m)*n(m))*gt(k));
%    e2 = e2 +...
%      (((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l)) * gn(l,k) + ...
%      (st*(n(i)*n(j) - dd(i,j))*dd(j,l) - ct*ee(i,j,l)*n(l))*gt(k)) *...
%      (((1-ct)*(dd(i,m)*n(k) + dd(k,m)*n(i)) - st*ee(i,k,m)) * gn(m,j) + ...
%      (st*(n(i)*n(k) - dd(i,k))*dd(k,m) - ct*ee(i,k,m)*n(m))*gt(j));
%  end; end; end; end; end;
%  return

  %%% variant 3 - to see all terms
  e1 = 0;
  e2 = 0;
  u=0;

  for k=1:3; for j=1:3; for i=1:3;
    u = u + (1-ct)^2 * n(i)*n(j)*gn(k,i)*gn(k,j);
  end; end; end

  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 - 2*st*(1-ct) * ee(k,j,l)*n(k)*gn(l,j)*gn(i,i);
    e2 = e2 - 2*st*(1-ct) * ee(k,j,l)*n(k)*gn(l,i)*gn(i,j);
  end; end; end; end

  for i=1:3; for j=1:3;
    u = u + st^2 * gn(i,j)*gn(i,j);
    e1 = e1 + (1-ct)^2 * gn(i,i)*gn(j,j);
    e2 = e2 - st^2 * gn(i,i)*gn(j,j);
    e1 = e1 - st^2 * gn(i,j)*gn(j,i);
    e2 = e2 + (1-ct)^2 * gn(i,j)*gn(j,i);
  end; end

  for i=1:3; for j=1:3;
    e1 = e1 - 2*st * n(i)*gn(j,i)*gt(j);
    e2 = e2 - 2*st * n(j)*gn(i,i)*gt(j);
  end; end

  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 + 2*ct*(1-ct) * ee(k,j,l)*n(k)*n(i)*gn(l,i)*gt(j);
    e2 = e2 - 2*st^2      * ee(k,j,l)*n(k)*n(i)*gn(l,i)*gt(j);
    e2 = e2 + 2*ct*(1-ct) * ee(k,i,l)*n(j)*n(k)*gn(l,i)*gt(j);
    e1 = e1 - 2*st^2      * ee(k,i,l)*n(k)*n(j)*gn(l,i)*gt(j);
  end; end; end; end

  for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 + 2*st^2 * ee(j,i,k)*gn(k,i)*gt(j);
    e2 = e2 - 2*st^2 * ee(j,i,k)*gn(k,i)*gt(j);
  end; end; end

  for j=1:3; for i=1:3;
     u=u+(dd(i,j)-n(i)*n(j)) * gt(i)*gt(j);
  end; end
  e1=e1+u; e2=e2+u;
  return

end

