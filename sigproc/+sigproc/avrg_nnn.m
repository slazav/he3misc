function [xs,ys,exs,eys] = avrg_N(x,y, groups)
  % average groups of points
  xs=[]; ys=[]; exs=[]; eys=[];
  n=1;
  for i=1:length(groups);
    ii=n:n+groups(i)-1;
    xs(end+1) = mean(x(ii));
    ys(end+1) = mean(y(ii));
    exs(end+1) = sqrt(sum((x(ii) - mean(x(ii))).^2) / length(ii));
    eys(end+1) = sqrt(sum((y(ii) - mean(y(ii))).^2) / length(ii));
    n=n+groups(i);
  end
  xs=xs';
  ys=ys';
  exs=exs';
  eys=eys';
end
