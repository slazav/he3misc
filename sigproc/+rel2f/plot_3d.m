function plot_3d(time, freq, amp, scale, fig_title)
% plot_3d(time, freq, amp, scale, title)
%
% plot 3d-like picture of a signal
% scale can be 'auto' or real number
%
% slazav, feb 2012

    %% plot 3d-like frequency picture 

    [p_f, p_a] = rel2f.find_peaks(freq,amp(:,1));
    [p2f, p2a] = rel2f.select_2peaks(p_f, p_a);

    f1 = min(p2f);
    f2 = max(p2f);
    f25 = round(f2/25)*25;

    if     strcmp(scale, 'auto'); s = 2/max(max(amp));
    elseif strcmp(scale, 'log');  s = 0.02/log(max(max(amp)));
    else   s=scale;
    end

    find_figure(['3d: ' fig_title]);
    hold off;
    for i=-2:3
      plot( [1,1]*i*25+f25, [time(1),time(length(time))], '-r'); hold on;
    end
    for i=1:length(time)
      if strcmp(scale, 'log')
        plot(freq, s*log(amp(:,i)) + time(i),'-b');
      else
        plot(freq, s*amp(:,i) + time(i),'-b');
      end
    end
    hold off;

    xlim( f2 + 300*[-1,1] );
    ylim( 15*[-1,1] );

    ht=title(fig_title);
    set(ht,'Interpreter','none');
end

