function [amp, f0 ] = integrate_peak(tx,xx, df, f0)
% [amp, f0 ] = integrate_peak(tx,xx, df, f0)
%
% Integrates signal inside the +/-df band around main peak
% in the f0+/-df region.
% Returns amplitude and peak position (w/o smoothing)

    N = length(tx);
    ddt = (tx(N)-tx(1))/(N-1);
    ddf = 1/ddt/N; % frequency step in fft

    X = fft(xx);

    % find maximum position f0 around old value
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    [a0m, f0i] = max(abs(X(ii)));
    f0=(ii(1) + f0i)*ddf;

    % filtering range for the fixed value of f0
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];

    amp = sqrt(sum(abs(X(ii) .* blackman(length(ii))').^2))/N;

%%%  more correct peak position, slow
%    af = sqrt( sum(ii.^2 .* abs(X(ii) .* blackman(length(ii))').^2))/N;
%    f0 = af/amp *ddf;

end
