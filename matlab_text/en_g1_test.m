function en_g1_test()
  figure;
  [a t] = meshgrid(0:10:180);
  o=ones(size(a));
  b = o*10;

  dx = o*0.1;
  dy = o*0.1;
  dz = o*0.1;

  dax = o*0.011;
  day = o*0.012;
  daz = o*0.013;

  dbx = o*0.014;
  dby = o*0.015;
  dbz = o*0.016;

  dtx = o*0.017;
  dty = o*0.018;
  dtz = o*0.019;

  E0 = arrayfun(@en_g1_0, a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz);
  E1 = arrayfun(@en_g1  , a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz);
  mesh(a,t, E1-E0);

  figure;
  [b t] = meshgrid(0:10:180);
  a = o*10;

  E0 = arrayfun(@en_g1_0, a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz);
  E1 = arrayfun(@en_g1  , a,b,t, dax,day,daz, dbx,dby,dbz, dtx,dty,dtz, dx,dy,dz);
  mesh(b,t, E1-E0);

end
