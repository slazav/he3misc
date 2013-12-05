function [f, a] = fft(t, x, minf, maxf)
%[f, a] = fft(t, x)
%
% Do FFT of x(t).

    if nargin < 3 minf=0; end
    if nargin < 4 maxf=2000; end

    N = length(t);
    a = 2*fft(x,N)/(N-1);

    f = linspace(0,1,N)*N/(t(end)-t(1));
    if maxf<=0; maxf=max(f)/2; end
    if minf<=0; minf=0; end
    ii = find(f >= minf & f < maxf);
    f=f(ii);
    a=a(ii);
end
