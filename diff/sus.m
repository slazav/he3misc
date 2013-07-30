function [ sus] = sus(P,TTc)
%chi, following D.Einzel JLTP 84 and
Fa0= F0a(P);


%calc sus
sus=sus0(P,TTc)/(1+Fa0*sus0(P,TTc));


end

