function en_gr_test1()
  [a b] = meshgrid(0:0.2:pi);
  o=ones(size(a));
  t = o*acos(-0.25);

  dx = o*1e-6;
  dy = o*1e-6;
  dz = o*1e-6;

  dax = o*0.11 * dx;
  day = o*0.12 * dx;
  daz = o*0;

  dbx = o*0.14 * dx;
  dby = o*0.15 * dx;
  dbz = o*0;

  dtx = o*0;
  dty = o*0;
  dtz = o*0;

  [E10 E20] = arrayfun(@en_gr, a,b,t, dax,a.*day,daz, dbx,dby,dbz, t.*dtx,dty,dtz, dx,dy,dz);
  [E1 E2]   = arrayfun(@en_gr1, a,b, dax,a.*day, dbx,dby, dx,dy);

  figure; hold on
  mesh(a,b, (E10-E1));
  mesh(a,b, E10);
  view(3);

  figure; hold on
  mesh(a,b, (E20-E2));
  mesh(a,b, E20);
  view(3);

end
