function yc = autocorr(X,Y,xc)

  Y = Y-mean(Y);

  ii=[];
  for i=2:length(X)
    if X(i)~=X(i-1); ii(end+1)=i; end
  end
  X=X(ii); Y=Y(ii);

  for i=1:length(xc)
    Y1 = interp1(X,Y, X+xc(i));
    ii = find(~isnan(Y1));
    yc(i) = sum(Y(ii) .* Y1(ii)) / length(ii);
  end
end
