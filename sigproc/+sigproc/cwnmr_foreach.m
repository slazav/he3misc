function cwnmr_foreach(func, time, fset, fmeas, Abs, Disp)
% run func for each sweep:
%   func(dir, time, fset, fmeas, Abs, Disp);
% use cwnmr_get() to get raw input data

  pdir=0; % previous sweep direction. 
  ppt=1;  % previous point of sweep direction change. 
  for i=2:length(time)
    % sweep direction in the point. 
    if     fset(i) > fset(i-1); d=1;
    elseif fset(i) < fset(i-1); d=-1;
    else d=0;
    end
    % direction change or last point.
    if d ~= pdir || i == length(time)
      func(pdir, time(ppt:i),...
        fset(ppt:i), fmeas(ppt:i), Abs(ppt:i), Disp(ppt:i));
      pdir=d;
      ppt=i-1;
    end
  end
end
