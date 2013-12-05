function [amp, f0 ] = integrate_peak_fft(X, ddf, df, f0)
% [amp, f0 ] = integrate_peak_fft(X, ddf, df, f0)
%
% the same as integrate_peak, but gets fft instead of raw signal

    % find maximum position f0 around old value
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    ii = ii(find(ii>0 & ii<length(X)));
    if (length(ii)<1); amp=0; return; end

    [a0m, f0i] = max(abs(X(ii)));
    f0=(ii(1) + f0i)*ddf;

    % filtering range for the fixed value of f0
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    ii = ii(find(ii>0 & ii<length(X)));

    amp = sqrt(sum(abs(X(ii) .* blackman(length(ii))').^2))/length(X);

%%%  more correct peak position, slow
%    af = sqrt( sum(ii.^2 .* abs(X(ii) .* blackman(length(ii))').^2))/N;
%    f0 = af/amp *ddf;

end
