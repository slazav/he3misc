#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  subplot(2,2,1); hold on;
  plot_tpdep(@he3_text_lg1, 'northeast');
  plot_tpdep(@he3_text_lg2);
  title('\lambda_{G1}, \lambda_{G2}, erg/cm vs T/T_c');
  text(0.50, 1e-11, '\lambda_{G2}');
  text(0.87, 2e-11, '\lambda_{G1}');
  ylim([0 6e-11])

  subplot(2,2,2); hold on;
  plot_tpdep(@he3_text_delta);
  title('\delta vs T/T_c');
  ylim([-0.35 0])

  subplot(2,2,3); hold on;
  plot_tpdep(@he3_text_cperp);
  title('c_\perp, cm/c vs T/T_c');
%  ylim([-0.35 0])

  subplot(2,2,4); hold on;
  plot_tpdep(@he3_text_c);
  title('c vs T/T_c');
  ylim([0 2.5e-10])

  print text_gr.eps -deps -color "-S1200,1000"
