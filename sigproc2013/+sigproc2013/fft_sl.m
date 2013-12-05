function [time, freq, amp, window, step] = fft_sl(tx,xx, window, step, minf, maxf)
% Do sliding FFT of xx(tx) with some window and step.
% Parallel calculations are used.
%
% -- slazav, feb 2012, dec 2013.

    if nargin < 3 || window<=0; window=floor(length(tx)/100); end
    if nargin < 4 || step<=0;   step=floor(window/10); end
    if nargin < 5 || minf<0; minf=10; end
    if nargin < 6 || maxf<0; maxf=3000; end

    fprintf('running fft (window: %d, step: %d)\n', window, step);
    wind = blackman(window)';
%    wind = ones(1,window);
    amp=[]; freq=[]; time=[];
    maxi = length(tx)-window;
    ii=(1:step:maxi);
    parfor j=1:length(ii)
      i=ii(j);
      [freq(:,j), amp(:,j)] =...
        sigproc2013.fft(tx(i:i+window-1), xx(i:i+window-1).*wind, minf, maxf);
      time(j)=tx(round(i+window/2)); %The recorded times are the middle of each window.
    end
    freq=freq(:,1);
end
