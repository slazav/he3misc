function [xs,ys] = avrg_dx(x,y, dx)
  % average data in y-dy/2 .. y+dy/2 range
  xs=[]; ys=[];
  for i=min(x):dx:max(x)
    ii=find(x > i-dx/2 & x <= i+dx/2);
    if length(ii)
      xs(end+1) = mean(x(ii));
      ys(end+1) = mean(y(ii));
    end
  end
end
