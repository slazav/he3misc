function rel = he3ppd_rel(Omega, Fork)
% PPD relaxation as a function of Omega and Fork1 width
% see /rota/Analysis/PS/osc2011/Relax/relax_temp

  rel = (0.0244244 + 0.114521 * Omega) + (0.000859956 + 0.000683031 * Omega) * Fork;

end