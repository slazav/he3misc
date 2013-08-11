function [ D ] = angle_integr_lowtemp( P,TTc,w0,ksi )
%Following following D.Einzel JLTP 84 AND Markelov & Mukharsky Physica B
%178

%PARAMETERS


Fa0 = he3_f0a(P);
gamma=he3_gyro;
H=w0/gamma;

w_exch=-Fa0*gamma*sus(P,TTc)*H; 


tauDp = tauDperp( TTc,P );

kB=8.6173*10^(-5); %eV/K
[T Tc]=TTc_to_T(P,TTc);

sgap=he3_trivgap(TTc,P)*kB*Tc;

Vf=he3_vf(P);
suspar=sus(P,TTc);

%INTEGRAL
scale=10^7;%exp(-sgap/(kB*Tc*TTc));  %quad uses absolute error tolerance, wipe out strongest temperature dep.

%testing
%  figure;hold on
%  the=0:pi/1000:pi;
%  for kk=1:1000
%      phi=2*pi*kk/1000;
%      int=integrand(the,phi,ksi,tauDp,w0,w_exch,sgap,kB*T,Vf,suspar,scale);
%      plot(the,int)
%  end
% drawnow;

Dt=2/(2*2*pi)*dblquad(@(the,phi) integrand(the,phi,ksi,tauDp,w0,w_exch,sgap,kB*T,Vf,suspar,scale),0,pi,0,pi*2,10^-9);

D=Dt*scale;  %quad uses absolute error tolerance, put back the main temperature dep.

end

function intgr= integrand(the,phi,ksi,tauDp,w0,w_exch,sgap,kBT,Vf,suspar,scale)

dOmega=sin(the);
kz=sin(the);  %z-axis chosen to kz-direction


s=(w0+w_exch)*tauDp/(1-1i*w0*tauDp);
integr_ksi = integrand_ksi(ksi,the,phi,kz,s,sgap,kBT);

intgr_help=dOmega.*tauDp./(1-1i*w0*tauDp).*kz.^2.*integr_ksi;
intgr=real(intgr_help)./(4*kBT)*Vf^2/suspar;
intgr=intgr./scale;  %quad uses absolute error tolerance, wipe out strongest temperature dep.
end

function integr_ksi = integrand_ksi(ksi,the,phi,kz,s,sgap,kBT)
    
    
    Ek=sqrt(ksi.^2+sgap^2);
    %phik=(cosh(Ek./(2*kBT))).^(-2);
    
    %Markelov & Mukharsky PRB 178
    %Sminus_sq=1+1/2*(1-kz^2)*(-1+ksi.^2/Ek.^2);
    %Splus_sq=ksi.^2./Ek.^2+kz.^2.*sgap.^2./Ek.^2;
    
    %Bunkov et al PRL 65
     kx=cos(the).*cos(phi);
     Sminus_sq=(ksi./Ek.*kx.^2+1-kx.^2).^2;
    Splus_sq=(ksi./Ek.*(1-kz.^2)+kz.^2).^2;
    
    %Einzel JLPT 84
    %Sminus_sq=(ksi./Ek.*kz.^2+1-kz.^2).^2;
    %Splus_sq=(ksi./Ek.*(1-kz.^2)+kz.^2).^2;
    
    help1=(Sminus_sq-1i*s*Splus_sq)./(1+s.^2.*Splus_sq);
    integr_ksi=help1;
end