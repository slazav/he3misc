function [time, freq, amp] = fft_sl(tx,xx, window, step, minf, maxf)
% [time, freq, amp] = fft_sl(tx,xx, window, step, [minf] [maxf])
%
% Do sliding FFT of xx(tx) with some window and step
%
% -- slazav, feb 2012.

    if nargin < 5; minf=10; end
    if nargin < 6; maxf=2000; end

    fprintf('running fft (window: %d, step: %d):   0 %%', window, step);
    j = 1;
    maxi = length(tx)-window;
    pr1=0;
    wind = blackman(window)';
%    wind = ones(1,window);
    time=[]; amp=[]; freq=[];
    for i=1:step:maxi
        pr2=int32(100*i/maxi);
        if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); pr1=pr2; end
        [freq, amp(:,j)] =...
          sigproc.fft(tx(i:i+window-1), xx(i:i+window-1).*wind, minf, maxf);
        time(j)=tx(round(i+window/2)); %The recorded times are the middle of each window.
        j=j+1;
    end
    fprintf('\b\b\b\b\bok   \n');
end
