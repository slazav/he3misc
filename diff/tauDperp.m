function [ tauDperp ] = tauDperp( TTc,P )
[ tau_aver,~] = tau(P,TTc);
%following D.Einzel JLTP 84, and D. Einzel JLTP 54

lambda=scattering_factors(P);

yos0=yosida(P,TTc,0);
yos2=yosida(P,TTc,2);

tauDperp=tau_aver/(1-lambda*(4/5*yos0+1/5*yos2)/yos0);
end


function [lambda]=scattering_factors(P)
%notation: Einzel JLTP54 (A1 and A0)
%following Einzel & Wölfle JLTP 32 page 34 and 27
%A1=At, A0=As  

[As0 As1 Aa0 Aa1] =landau_pars(P);

S1=As1-3*Aa1;
S0=As0-3*Aa0;

T1=As1+Aa1;
T0=As0+Aa0;

%averages over Abrikosov angles, see Collision_integrals.nb
Waver=pi/60*(30*S0^2-20*S0*S1+14*S1^2+45*T0^2-30*T0*T1+21*T1^2);

W=1/420*pi*(-70*S0^2-54*S1^2+175*T0^2+28*S0*(7*S1+10*T0-6*T1)-42*T0*T1+71*T1^2+8*S1*(-21*T0+19*T1));


lambda=W/Waver;
end


function [As0 As1 Aa0 Aa1] =landau_pars(P)
%Dobbs page 52 table 3.2 (It is not clear what parameters Einzel used in
%his derivation of the spin diffusion coefficient. These parameters
%slightly disagree with table 1 in Einzel & Wölfle JLTP 32)

Ptable=[0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3 3.3 3.439]*10; %bar

As0_tab=[0.903 0.941 0.957 0.967 0.972 0.977 0.980 0.982 0.984 0.986 0.987 0.988 0.989];
As1_tab=[1.927 2.052 2.139 2.204 2.256 2.300 2.338 2.372 2.403 2.430 2.455 2.477 2.487];
Aa0_tab=-[2.311 2.623 2.759 2.861 2.937 3.016 3.065 3.098 3.115 3.115 3.098 3.098 3.082];
Aa1_tab=-[0.685 0.955 1.055 1.134 1.241 1.286 1.350 1.395 1.424 1.440 1.477 1.447 1.431];

As0=interp1(Ptable,As0_tab,P);%,'spline');
As1=interp1(Ptable,As1_tab,P);%,'spline');
Aa0=interp1(Ptable,Aa0_tab,P);%,'spline');
Aa1=interp1(Ptable,Aa1_tab,P);%,'spline');
end