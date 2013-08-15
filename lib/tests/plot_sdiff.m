function plot_sdiff()
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0:0.01:1;
  p=30;

  plot(ttc, he3_sdiff_nh(ttc, p), 'r-');
  plot(ttc, he3_sdiff_hperp(ttc, p), 'b-');
  plot(ttc, he3_sdiff_hpar(ttc, p), 'g-');
  plot(ttc, he3_sdiff(ttc, p, 1e6), 'c-', 'linewidth', 2);
  plot(ttc, he3_sdiff(ttc, p, 1e4), 'b-', 'linewidth', 2);
  plot(ttc, he3_sdiff(ttc, p, 0), 'm--');
  xlim([0 1])
  ylim([0 0.4])

  % improved Samuli's code
  %addpath ~/he3lib/diff
  %s1d=diff_coeff(p, ttc, 1e4);
  %s2d=diff_coeff(p, ttc, 1e6);
%  s2=load('sdiff_sam2');
%  plot(s2.ttc, s2.s1d/dn, 'm-');
%  plot(s2.ttc, s2.s2d/dn, 'm-');

  legend(
   'Normal hydrodynamic',...
   'Superfluid hydrodynamic D_{perp}',...
   'Superfluid hydrodynamic D_{par}',...
   'Superfluid D_{perp} 1000 kHz',...
   'Superfluid D_{perp} 0 Hz',...
  '');

%   'Superfluid D_{perp} 10 kHz',...

%  print -deps -color plot_sdiff.eps
end