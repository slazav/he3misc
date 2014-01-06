function e = en_g2_0(a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz)

  r = abt2r(a, b, t);
  dr(:,:,1) = (abt2r(a+dax, b+dbx, t+dtx) - r)/dx;
  dr(:,:,2) = (abt2r(a+day, b+dby, t+dty) - r)/dy;
  dr(:,:,3) = (abt2r(a+daz, b+dbz, t+dtz) - r)/dz;

  e=0;
  for k=1:3; for j=1:3; for i=1:3;
    e = e + dr(i,j,k) * dr(i,k,j);
  end; end; end
end

