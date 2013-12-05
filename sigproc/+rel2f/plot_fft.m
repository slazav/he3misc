function plot_fft(tx,xx, t0, window, fig_title)
% plot_fft(tx(i),xx(i), t0, window, fig_title)
%
% plot FFT at some fixed value
% scale can be 'auto' or real number
%
% slazav, feb 2012

    dt=(tx(length(tx))-tx(1))/(length(tx)-1);
    ii1=round((t0-tx(1))/dt)+1;
    if ii1<1 ii1=1; end
    ii2=ii1 + round(window/dt);
    if ii2>length(tx) ii2=length(tx); end

    [freq1, ar, ai] = rel2f.fft_ri(tx(ii1:ii2), xx(ii1:ii2));
    aa=sqrt(ai.^2+ar.^2);
    [~, mi]=max(aa);

    [freq2, ar2, ai2] = ...
      rel2f.fft_prec(tx(ii1:ii2), xx(ii1:ii2), freq1(mi), 60, 0.05);
    aa2=sqrt(ai2.^2+ar2.^2);

    find_figure(['fft: ' fig_title]); hold off;
    plot(freq1, aa,'-or');  hold on;
    plot(freq1, ar-0.3,'-ob');
    plot(freq1, ai-0.8,'-og');

    plot(freq2, aa2,'-m');
    plot(freq2, ar2-0.3,'-k');
    plot(freq2, ai2-0.8,'-k');

    xlim(freq1(mi) + 60*[-1,1]);

    ht=title(fig_title);
    set(ht,'Interpreter','none');
end