function [ Fa0] = F0a(P)
%parameter F_0^a, following D.Einzel JLTP 84

Fa0_vs_P=-[0.7 0.74 0.76 0.75];
Fa0P=[0 10 20 30];
Fa0=interp1(Fa0P,Fa0_vs_P,P,'spline');

end

