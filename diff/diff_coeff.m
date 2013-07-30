function [ D ] = diff_coeff( P,TTc,w0 )
%Following following D.Einzel JLTP 84 AND Markelov & Mukharsky Physica B
%178

%PARAMETERS


Fa0 = F0a(P);
gamma=2.03780000*10^8;%in Teslas
gamma=gamma*10^-4;% to CGS
H=w0/gamma;

w_exch=-Fa0*gamma*sus(P,TTc)*H; 


tauDp = tauDperp( TTc,P );

kB=8.6173*10^(-5); %eV/K
[T Tc]=TTc_to_T(P,TTc);

addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Spin_diffusion/gap')
sgap=trivial(TTc,P)*kB*Tc;

Vf=vff(P)*10^3;%cm/s
suspar=sus(P,TTc);

%INTEGRAL
%find relevant energy scale
phimax=(cosh(sgap/(2*kB*T)))^(-2);
steps=4;
phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);

while phimax<2*10^6*phi_steps_of_gap
    steps=steps*1.2;
 phi_steps_of_gap=(cosh(sgap*steps/(2*kB*T)))^(-2);
end

scale=exp(-sgap/(kB*Tc*TTc));  %quad uses absolute error tolerance, wipe out strongest temperature dep.


Dt=2/(2*2*pi)*triplequad(@(ksi,the,phi) integrand(ksi,the,phi,tauDp,w0,w_exch,sgap,kB*T,Vf,suspar,scale),0,sgap*steps,0,pi,0,pi*2,10^-7);

%testing
% ksit=0:sgap*steps/100:sgap*steps;
% figure;hold on
% for kk=1:100
%     thet=pi*kk/100;
%     int=integrand(ksit,thet,tauDp,w0,w_exch,sgap,kB*T);
%     [int(1) int(end)]
%     plot(ksit,int)
% end

D=Dt*scale;  %quad uses absolute error tolerance, put back the main temperature dep.

end

function intgr= integrand(ksi,the,phi,tauDp,w0,w_exch,sgap,kBT,Vf,suspar,scale)

dOmega=sin(the);
kz=sin(the);  %z-axis chosen to kz-direction


s=(w0+w_exch)*tauDp/(1-1i*w0*tauDp);
integr_ksi = integrand_ksi(ksi,the,phi,kz,s,sgap,kBT);

intgr_help=dOmega*tauDp/(1-1i*w0*tauDp)*kz^2*integr_ksi;
intgr=real(intgr_help)/(4*kBT)*Vf^2/suspar;
intgr=intgr/scale;  %quad uses absolute error tolerance, wipe out strongest temperature dep.
end

function integr_ksi = integrand_ksi(ksi,the,phi,kz,s,sgap,kBT)
    
    
    Ek=sqrt(ksi.^2+sgap^2);
    phik=(cosh(Ek./(2*kBT))).^(-2);
    
    %Markelov & Mukharsky PRB 178
    %Sminus_sq=1+1/2*(1-kz^2)*(-1+ksi.^2/Ek.^2);
    %Splus_sq=ksi.^2./Ek.^2+kz.^2.*sgap.^2./Ek.^2;
    
    %Bunkov et al PRL 65
     kx=cos(the)*cos(phi);
     Sminus_sq=(ksi./Ek.*kx.^2+1-kx.^2).^2;
    Splus_sq=(ksi./Ek.*(1-kz.^2)+kz.^2).^2;
    
    %Einzel JLPT 84
    %Sminus_sq=(ksi./Ek.*kz.^2+1-kz.^2).^2;
    %Splus_sq=(ksi./Ek.*(1-kz.^2)+kz.^2).^2;
    
    help1=(Sminus_sq-1i*s*Splus_sq)./(1+s.^2.*Splus_sq);
    integr_ksi=phik.*help1;
end