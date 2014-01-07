function en_gr_test()
  [a t] = meshgrid(0:0.2:pi);
  o=ones(size(a));
  b = o*0.1;

  dx = o*1e-6;
  dy = o*1e-6;
  dz = o*1e-6;

  dax = o*0.11 * dx;
  day = o*0.12 * dx;
  daz = o*0.13 * dx;

  dbx = o*0.14 * dx;
  dby = o*0.15 * dx;
  dbz = o*0.16 * dx;

  dtx = o*0.17 * dx;
  dty = o*0.18 * dx;
  dtz = o*0.19 * dx;

  [E10 E20] = arrayfun(@en_gr_0, a,b,t, dax,a.*day,daz, dbx,dby,dbz, t.*dtx,dty,dtz, dx,dy,dz);
  [E1 E2]   = arrayfun(@en_gr, a,b,t, dax,a.*day,daz, dbx,dby,dbz, t.*dtx,dty,dtz, dx,dy,dz);

  figure; hold on
  mesh(a,t, (E10-E1));
  mesh(a,t, E10);
  view(3);

  figure; hold on
  mesh(a,t, (E20-E2));
  mesh(a,t, E20);
  view(3);

  [b t] = meshgrid(0:0.2:pi);
  a = o*0.1;

  figure; hold on
  [E10 E20] = arrayfun(@en_gr_0, a,b,t, dax,day,daz, b.*dbx,dby,dbz, dtx,dty,t.*dtz, dx,dy,dz);
  [E1 E2]   = arrayfun(@en_gr, a,b,t, dax,day,daz, b.*dbx,dby,dbz, dtx,dty,t.*dtz, dx,dy,dz);
  mesh(b,t, (E10-E1));
  mesh(b,t, E10);
  view(3);

  figure; hold on
  mesh(b,t, (E20-E2));
  mesh(b,t, E20);
  view(3);

end
