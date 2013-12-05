function [time, freq] = findfreq(tx,xx)

  % zero crossing times
  tc = rel2f.find_zeros(tx,xx);

  if 0
    find_figure('findfreq'); clf; hold on;
    plot(tx,xx,'o-k','MarkerSize',3);

    plot(tc,zeros(size(tc)),'sm');
    plot([tx(1) tx(end)],[0 0],'-','Color',[0.6 0.6 0.6]);
    xlabel('time, s');
    ylabel('filtered signal');
    set(gcf,'renderer','zbuffer');
  end

  time = (tc(1:end-1)+tc(2:end))/2;
  freq = 0.5 ./ diff(tc);
end
