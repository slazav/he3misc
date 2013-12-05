function [amp, fre, dis, f0 ] = integrate_peakf_fft(X, ddf, df, f0)
% fre = integrate_peakf_fft(X, ddf, df, f0)
%
% 1st moment of the peak, gets fft instead of raw signal

    % find maximum position f0 around old value
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    ii = ii(find(ii>0 & ii<length(X)));
    [a0m, f0i] = max(abs(X(ii)));
    f0=(ii(1) + f0i)*ddf;

    % filtering range for the fixed value of f0
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    ii = ii(find(ii>0 & ii<length(X)));

    power=abs(X(ii)).^2 .* blackman(length(ii))';

    m0 = sum(power);
    m1 = sum(ii .* power);
    m2 = sum((ii-m1/m0).^2 .* power);

    amp = sqrt(m0)/length(X);
    fre = m1/m0 *ddf;
    dis = sqrt(m2/m0) *ddf;

%    amp = sqrt(sum(abs(X(ii) .* blackman(length(ii))').^2))/length(X);



end
