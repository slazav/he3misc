function [xs,ys] = avrg_dx(x,y, dx)
  % average data in y-dy/2 .. y+dy/2 range
  xs=[]; ys=[];
  i=1;
  while i<=length(x)
    ii=find(x-x(i) >= 0 & x-x(i) < dx);
    if length(ii)
      xs(end+1) = mean(x(ii));
      ys(end+1) = mean(y(ii));
    end
    i=ii(end)+1;
  end
end
