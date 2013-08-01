function texture_05bar

%  [textur,spec] = texture(ndisc,textpar,nsp,specpar,initype)
%
%  1st input argument is the number of discretization intervals
%  for texture, or empty (then do not optimize the texture,
%  use initial value)
%  
%  2nd input argument is texture calculation parameters -
%  vector with elements:
%  1 - temperature / Tc
%  2 - pressure, bar
%  3 - larmor frequency, kHz
%  4 - cylinder radius, cm
%  5 - rotation velocity, rad/s
%  6 - omega_v, rad/s
%  7 - lambda/omega (s/rad) if >=0
%      use calculated lambda/omega if == -1
%  8 - lambga_HV (kg/(m^3 T^2)) if >= 0
%      use calculated lambga_HV if == -1
%  
%  3rd input argument is the number of intervals in the output spectrum,
%  or empty if spectrum is not needed
%  
%  4th input argument is spectrum calculation parameters -
%  vector with elements :
%  1 - half-width of NMR line
%  2 - margin for automatic region determination
%  
%  5th input argument is initial condition
%  1 - texture w/o 90 deg peak
%  2 - texture with 90 deg peak
%  n by 3 array - 1st column ignored , 2nd - alpha, 3rd - beta
%  
%  1st output argument
%  texture - n by 3 array - 1st column r, 2nd - alpha, 3rd - beta
%  
%  2nd output argument
%  spectrum - m by 2 array - 1st column freq, 2nd column absorption
