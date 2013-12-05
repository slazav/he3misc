function [p_f, p_a] = find_peaks_sl(freq, amp)
% [p_f(k,t), p_a(k,t)] = find_peaks_sl(freq(n), amp(n,t))
%
% Find peaks for all time values.
% Use output of fft_sl() as input.
%
% -- slazav, feb 2012.

    fprintf('finding peaks:   0 %%');
    maxk=length(amp(1,:));
    for k=1:maxk
      pr1=int32(100*(k-1)/maxk);
      pr2=int32(100*k/maxk);
      if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); end
      mm=1;
      [max_f, max_a] = rel2f.find_peaks(freq, amp(:,k));
      p_f(1:length(max_f), k) = max_f;
      p_a(1:length(max_a), k) = max_a;
    end
    fprintf('\b\b\b\b\bok   \n');
end
