function show_fft(dstr, xfile, t1, t2)
% show_fft(dstr, xfile, t1, t2)

%    if nargin < 2; [dstr, xfile] = rel2f.last_sig(); end

    fig_title = [dstr ' ' xfile];
    [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile);
    fprintf('dt_osc: %g, points: %d\n', dt_osc, length(tx));

    ii = find(tx>str2num(t1) & tx<=str2num(t2));

    [fre,amp] = sigproc.fft(tx(ii), xx(ii),2,5000);
    find_figure(fig_title); hold off;
    plot(fre,abs(amp),'-r'); hold on;
%    xlim([1,300]);
end