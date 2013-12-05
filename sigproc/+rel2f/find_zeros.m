function tz = find_zeros(tx,xx)
% find zero crossings

  xxp = xx > 0;
  xxm = ~xxp;
  jc = find((xxp(1:end-1) & xxm(2:end)) | (xxm(1:end-1) & xxp(2:end)));

  t1 = tx(jc); t2 = tx(jc+1);
  x1 = xx(jc); x2 = xx(jc+1);

  % crossing times
  tz = t1 - x1.*(t2-t1)./(x2-x1);
end
