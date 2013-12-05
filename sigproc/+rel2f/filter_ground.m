function [ xxf, f0 ] = filter_ground(tx,xx, df, f0)
% [ xxf, f0 ] = filter_ground(tx,xx, df, f0)
%
% Remove frequencies outside the +/-df band around ground peak.
% Ground peak is a one of two main peaks with higher frequency.
% get approx value for f0, returns exact value

    plot=0;

    N = length(tx);
    ddt = (tx(N)-tx(1))/(N-1);
    ddf = 1/ddt/N; % frequency step in fft

    X = fft(xx);

    if plot
      ff=gcf;
      find_figure('freq filter'); clf; hold on;
      plot ([1:N]*ddf, abs(X), '-b');
      plot ([1,1]*f0-df, [0 max(abs(X))], '-k');
      plot ([1,1]*f0+df, [0 max(abs(X))], '-k');
    end

    % find maximum position f0 around old value
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    [~, f0i] = max(X(ii));
    f0=(ii(1) + f0i)*ddf;


    % apply filter
    ii = [round((f0-df)/ddf):round((f0+df)/ddf)];
    X1 = X(ii);
    X = zeros(1,N);
    X(ii) = X1 .* blackman(length(X1))';

    if plot
      plot (ddf*ii, abs(X(ii)), '-r');
      ylim([0 2*max(abs(X(ii)))]);
      figure(ff);
    end

    % reconstruct signal
    xxf = ifft(X);
end
