function [e1 e2] = en_gr_0(a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz)

  r = abt2r(a, b, t);
  gr(:,:,1) = (abt2r(a+dax, b+dbx, t+dtx) - r)/dx;
  gr(:,:,2) = (abt2r(a+day, b+dby, t+dty) - r)/dy;
  gr(:,:,3) = (abt2r(a+daz, b+dbz, t+dtz) - r)/dz;

  e1=0;
  for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 + gr(i,j,j) * gr(i,k,k);
  end; end; end

  e2=0;
  for k=1:3; for j=1:3; for i=1:3;
    e2 = e2 + gr(i,j,k) * gr(i,k,j);
  end; end; end

end

