function plot_amps(time, pf, pa, fig_title)
% plot_amps(time(t), pf(k,t), pa(k,t), fig_title)
%
% plot amplitudes of main peaks
%
% slazav, feb 2012

    [p2f, p2a] = rel2f.select_2peaks_sl(pf, pa);
    [psf, psa] = rel2f.select_side_peaks_sl(pf, pa, p2f(1,:), 25, 2);

    %% plot amplitude picture 
    find_figure(['amp: ' fig_title]);
    hold off;
    plot([0,0], [0,max(max(p2a))], '-k');
    hold on
    plot(time, p2a(1,:),'-r');
%    plot(time, p2a(2,:),'-c');
    plot(time, psa,'-g');
    hold off;

    ht=title(fig_title);
    set(ht,'Interpreter','none');
end