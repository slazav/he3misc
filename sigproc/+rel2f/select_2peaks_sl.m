function [p2_f, p2_a] =select_2peaks_sl(p_f, p_a)
% [p2_f(2,t), p2_a(2,t)] =select_2peaks_sl(p_f(t), p_a(t))
%
% Select two main peaks for all time values
% Use output from find_peaks_sl() as input.
%
% -- slazav, feb 2012.

    maxk=length(p_a(1,:));
    for k=1:maxk
      %just two largest peaks
      [p2_f(1:2,k), p2_a(1:2,k)] = rel2f.select_2peaks(p_f(:,k), p_a(:,k));
      % process noisy tails
      if (k<2) continue; end
      for n=1:2
        % if peak jumps, try to find nearest peak
        [ko, ~] = max([k-20,1]);
        fpr = predict_next(p2_f(n,ko:k-1));

        if abs(p2_f(n,k)-fpr)>10
          [~, near_idx] = min(abs(p_f(:,k)-fpr));
          p2_f(n,k)=p_f(near_idx,k);
          p2_a(n,k)=p_a(near_idx,k);
        end
%        % if it doesnt work, cut the signal off
%        if (abs(p2_f(n,k)-fpr)>10) || (p2_a(n,k-1)==0)
%          % cut to previous value
%          p2_f(n,k)=p2_f(n,k-1);
%          p2_a(n,k)=0;
%        end
      end
    end
end

% Liner prediction of the next value of the array arr.
% Used in select_2peaks_t()
function res = predict_next(arr)
  n = length(arr);
  sx  = n*(n+1)/2; % sum([1:n])
  sxx = sx*(2*n+1)/3; % sum([1:n].^2)
  sf  = sum(arr);
  sfx = sum(arr.*[1:n]);
  B = (sxx*sf-sx*sfx)/(sxx*n-sx*sx);
  A = (sf - B*n)/sx;
  res=A*(n+1) + B;
end
