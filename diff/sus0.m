function sus0 = sus0(P,ttc)
%chi0, following D.Einzel JLTP 84
% see also WV book ch.10 p.449
gap = he3_trivgap(ttc,P);
sus0 = 2/3 + 1/3 * he3_yosida(ttc,gap,0);
end

