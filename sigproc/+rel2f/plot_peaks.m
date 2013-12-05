function plot_peaks(time, pf, pa, fig_title)
% plot_peaks(time(t), pf(k,t), pa(k,t), fig_title)
%
% plot frequency picture with all maximums and detected peaks 
%
% slazav, feb 2012

    [p2f, p2a] = rel2f.select_2peaks_sl(pf, pa);
    [psf, psa] = rel2f.select_side_peaks_sl(pf, pa, p2f(1,:), 25, 2);

    f1 = min(min(p2f));
    f2 = max(max(p2f));
    f25 = round(f2/25)*25;

    find_figure(['peaks: ' fig_title]);
    hold off;
    plot(pf,time, '.r', 'MarkerSize', 0.5); hold on;
%    plot(time, p_f+25,'.m', 'MarkerSize', 0.5); 
    plot(p2f,time, '-b');
%    plot(time, psf,'-g');
    for i=-2:2
      plot(  [1,1]*i*25+f25, [time(1),time(length(time))], '-m');
    end
    hold off;

    xlim(f2 + 100*[-1,1] );

    ht=title(fig_title);
    set(ht,'Interpreter','none');
end