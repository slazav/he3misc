function plot_diff()
  find_figure('spin_diff'); clf; hold on;
  p=[0 34];
  for i=1:length(p)
    ttc=0:0.05:1;
    w=9200;
    plot(ttc, diff_coeff(p(i),ttc,w),'r-');
  end
end