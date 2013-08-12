function [ D ] = diff_coeff( P,TTc,w0 )
%Following following D.Einzel JLTP 84 AND Markelov & Mukharsky Physica B 178

  for i=1:length(TTc)
    Fa0 = he3_f0a(P);
    w_exch=-Fa0*w0*sus(P,TTc(i)); 
    tauDp = tauDperp( TTc(i),P );

    sgap=he3_trivgap(TTc(i),P);
    Vf=he3_vf(P); %cm/s
    suspar=sus(P,TTc(i));

    D(i) = dblquad(...
      @(x,the) integrand(x,the,tauDp,w0,w_exch,sgap,TTc(i),Vf,suspar),...
      0,1,0,pi,10^-7);
  end
end

function intgr= integrand(x,the,tauDp,w0,w_exch,sgap,TTc,Vf,suspar)

  ksi = atanh(x) *(2*TTc);

  kz = sin(the);

  t = tauDp/(1-1i*w0*tauDp);
  s = (w0+w_exch) * t;

  Ek = sqrt(ksi.^2+sgap^2);
  phik = (cosh(Ek./(2*TTc))).^(-2);
  u = ksi./Ek;

  Splus_sq = (u.*(1-kz.^2)+kz.^2).^2;

  help1=(3/8*(u-1).^2 * cos(the).^4 + (u-1)*cos(the).^2 +...
     1 - 1i*s*Splus_sq)./(1+s.^2.*Splus_sq);

  intgr_help = kz.^3 .* phik .* real(t.*help1);

  intgr=intgr_help/2 * Vf^2/suspar .* cosh(ksi/(2*TTc)).^2;
end


%>> diff_coeff(10, 0.4, 600000)
%ans = 0.522625
%ans = 0.522628 after setting correct range
