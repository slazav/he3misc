function [absamp] = fft_sl_remconst(amp, sm)
% Remove costant lines ftom signal.
% Returns abs value!

  for j=1:length(amp(:,1))
% the best method is to use 'rloess' smooth, but it is too slow...
%    absamp(j,:)=abs(abs(amp(j,:))-min(smooth(abs(amp(j,:)), 0.6, 'rloess')));
    absamp(j,:)=abs(abs(amp(j,:))-min(smooth(abs(amp(j,:)), sm)));
   end
end
