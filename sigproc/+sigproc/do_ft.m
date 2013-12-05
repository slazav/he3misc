function yc = do_ft(X1, Y1, xc)
  for i=1:length(xc)
    t = xc(i);
    fs=sin(2.0*pi/t * X1);
    fc=cos(2.0*pi/t * X1);
    yc(i) = hypot(sum(fs.*Y1), sum(fc.*Y1));
  end
end
