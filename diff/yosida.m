function [ yos,steps ] = yosida(P,TTc,n)
  %IMPLEMENTATION OF YOSIDA FUNCTIONS, FOLLOWING D.Einzel JLTP 84
  sgap=he3_trivgap(TTc,P);
  integ=2*quad(@(x) integrand(x,sgap,TTc,n),0,1,10^-8);
  yos=integ/(4*TTc);
end

function int = integrand(x,sgap,TTc,n)
    ksi = atanh(x);
    Ek=sqrt(ksi.^2+sgap.^2);
    phi=(cosh(Ek./(2*TTc))).^(-2);
    int=abs(ksi./Ek).^n .* phi .* cosh(ksi).^2;
end


