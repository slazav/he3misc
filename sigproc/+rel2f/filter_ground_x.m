function [ xxf, f0, f0ran ] = filter_ground_x(tx,xx, df, f0)
% [ xxf, f0 ] = filter_ground(tx,xx, df, f0)
%
% Remove frequencies outside the +/-df band around ground peak.
% Ground peak is a one of two main peaks with higher frequency.
% get approx value for f0, returns exact value

    N = length(tx);
    ddt = (tx(N)-tx(1))/(N-1);
    ddf = 1/ddt/N; % frequency step in fft
    X = fft(xx);
    if length(df) == 1
      df = [df df];
    end

    % find maximum position f0 around old value
    ii = [round((f0-df(1))/ddf):round((f0+df(2))/ddf)];
    [~, f0i] = max(X(ii));
    f0i=ii(1) + f0i -1;
    f0 = f0i*ddf;

    f0ran = [f0-df(1) f0+df(2)]; 

    % apply filter
    ii = [round((f0-df(1))/ddf):round((f0+df(2))/ddf)];
    X1 = zeros(1,N);

    X1(ii) = X(ii) .* blackman(length(ii))';

      find_figure('filter ground state'); clf
      plot (ddf*(1:N), abs(X), '-g'); hold on;
      plot (ddf*(1:N), abs(X1), '-r');
      plot (f0,  abs(X(f0i)), 'ob');
      xlim([f0-df(1)*1.2,f0+df(2)*1.2]);
      xlabel('signal frequency, Hz');
      ylabel('signal spectrum');
      legend('original','filtered','peak','Location','NorthEast');
      set(gcf,'renderer','zbuffer');

    % reconstruct signal
    xxf = ifft(X1);
end
