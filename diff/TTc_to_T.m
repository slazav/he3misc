function [ T,Tc] = TTc_to_T(P,TTc )
%T_to_Tc vs pressure, from thuneberg excel/Greywall86

Tc0=[2.438 2.425 2.411 2.395 2.378 2.360 2.339 2.317 2.293 2.267 2.239 2.209 2.177 2.143 2.106 ...
    2.067 2.026 1.981 1.934 1.883 1.828 1.769 1.705 1.636 1.560 1.478 1.388 1.290 1.181 1.061 0.997]*10^-3;   %Kelvins, Thuneberg spreadsheet/Greywall86
PTc=[30:-1:1 0.5]; %bar

Tc0=Tc0(end:-1:1);PTc=PTc(end:-1:1);

Tc=interp1(PTc,Tc0,P,'spline');

T=TTc*Tc;
end

