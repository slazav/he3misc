function h = show_fft3di(dstr, xfile)
%output is the handle for 3d-like image of the data

    cut_zero = 10; %Hz
    minf = 825000;
    maxf = 830000;

    if nargin < 2; [dstr, xfile] = sigproc.osc_last(); end

    fig_title = [dstr ' ' xfile];
    [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile);

    window=floor(length(tx)/100);
    step=floor(window/10);

%    window=round(0.8/dt_osc);
%    step=round(0.1/dt_osc);

    [time, freq, amp] = sigproc.fft_sl(tx, xx, window, step, minf, maxf);
    xx = []; tx = [];
    amp=abs(amp);
    if cut_zero>0; amp(find(freq<cut_zero),:) = 0; end

%amp=sigproc.fft_sl_remconst(amp);

    [amp] = sigproc.fft_sl_remconst(amp, 20);
%    [amp(:,find(time<=0))] = sigproc.fft_sl_remconst(amp(:,find(time<=0)), 20);

    find_figure(['3d: ' fig_title]); clf; hold on; title(fig_title);
    h = sigproc.plot_3di(time, freq, amp, 'sqrt');
%    h = sigproc.plot_3di(time, freq, amp, '');

%    unix(['mkdir -p -m 775 -- ', dstr, '/3d']);
%    fbase=[dstr, '/3d/' xfile];
%    print('-dpng', '-r150', [fbase, '.3d.png']);
%    hgsave([fbase, '.3d.fig']);
end

