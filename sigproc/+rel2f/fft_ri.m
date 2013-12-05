function [f, ar, ai] = fft_ri(t, x)
%[f, a, ph] = fft(t, x)
%
% Do FFT of x(t).

    L = length(t);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    X = fft(x,NFFT)/L;
    samplerate = 1/(t(2)-t(1));
    f = samplerate/2*linspace(0,1,NFFT/2+1);
    ar = real(X(1:NFFT/2+1));
    ai = imag(X(1:NFFT/2+1));
    j = find(f<=2000); %Some high enough frequency to limit the number of points
    f=f(j);
    ar=ar(j);
    ai=ai(j);
end
