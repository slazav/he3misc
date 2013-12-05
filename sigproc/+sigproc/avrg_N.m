function [xs,ys] = avrg_N(x,y, N)
  % average every N points
  xs=[]; ys=[];
  for i=1:N:length(x)-N+1;
    xs(end+1) = mean(x(i:i+N-1));
    ys(end+1) = mean(y(i:i++N-1));
  end
  xs=xs';
  ys=ys';
end
