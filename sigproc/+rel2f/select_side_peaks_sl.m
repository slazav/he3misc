function [ps_f, ps_a] = select_side_peaks_sl(p_f, p_a, f0, df, ddf)
% [ps_f(2,t), ps_a(2,t)] =
%        select_side_peaks_sl(p_f(n,t), p_a(n,t), f0(f), df,ddf)
%
% Select side peaks at f0+df +/- ddf and f0-df +/- ddf for
% all time values.
%
% -- slazav, feb 2012

    for k=1:length(p_f(1,:))
      [ps_f(1:2,k), ps_a(1:2,k)] = ...
        rel2f.select_side_peaks(p_f(:,k), p_a(:,k), f0(k), df, ddf);
%      ps_f(1:2,k)=f(1:2);
%      ps_a(1:2,k)=a(1:2);
    end
end
