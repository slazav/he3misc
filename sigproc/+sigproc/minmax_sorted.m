%% Sort peaks in a strange way:
%%  by minimal height difference with the lowest minimum
%%  on the way to the higher peak.
%% Code, which is working, but not used.
%% I don't want to lost it -- slazav

function [imin,imax, vmin,vmax] = minmax_sorted(a)
% find local min and max values, sort it 
  N=length(a);
  % finding local min and max indexes 
  imax=[]; imin=[];
  for i = 2:N-1.
    if a(i)>a(i-1) && a(i)>=a(i+1) imax(end+1) = i; end
    if a(i)<=a(i-1) && a(i)<a(i+1) imin(end+1) = i; end
  end
  % find weights for sorting: 
  % for max it is difference between its value and min in the way 
  % to higher or equal max 
  for i = 1:length(imax)
    vl = inf;
    for il=i-1:-1:1
      if a(imax(il)) >= a(imax(i)) % found larger or equal max on the left 
        iim = find(imin > imax(il) & imin < imax(i)); % minimums between these two max 
        vl = a(imax(i)) - min(a(imin(iim)));
        break;
      end
    end
    vr = inf;
    for ir=i+1:length(imax)
      if a(imax(ir)) >= a(imax(i)) % found larger or equal max on the right 
        iim = find(imin > imax(i) & imin < imax(ir)); % minimums between these two max 
        vr = a(imax(i)) - min(a(imin(iim)));
        break;
      end
    end
    vmax(i) = min(vl,vr);
    if isinf(vmax(i)); vmax(i) = a(imax(i)) - min(a(imin)); end
  end

  % the same for min 
  for i = 1:length(imin)
    vl = inf;
    for il=i-1:-1:1
      if a(imin(il)) <= a(imin(i)) % found smaller or equal min on the left 
        iim = find(imax > imin(il) & imax < imin(i)); % maxs between these two min 
        vl = max(a(imax(iim))) - a(imin(i));
        break;
      end
    end
    vr = inf;
    for ir=i+1:length(imin)
      if a(imin(ir)) <= a(imin(i)) % found smaller or equal min on the right 
        iim = find(imax > imin(i) & imax < imin(ir)); % maxs between these two min 
        vr = max(a(imax(iim))) - a(imin(i));
        break;
      end
    end
    vmin(i) = min(vl,vr);
    if isinf(vmin(i)); vmin(i) = max(a(imax)) - a(imin(i)); end
  end

  [vmin ii] = sort(vmin,'descend'); imin = imin(ii);
  [vmax ii] = sort(vmax,'descend'); imax = imax(ii);
end
