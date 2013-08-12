function [ yos,steps ] = yosida(P,TTc,n)
  %IMPLEMENTATION OF YOSIDA FUNCTIONS, FOLLOWING D.Einzel JLTP 84
  gap=he3_trivgap(TTc,P);
  for i=1:length(TTc)
    integ(i) = 2*quad(@(x) integrand(x,gap(i),TTc(i),n),0,1,10^-8);
  end
  yos=integ./(4*TTc);
end

function int = integrand(x,sgap,TTc,n)
    ksi = atanh(x) *(2*TTc);
    Ek=sqrt(ksi.^2+sgap.^2);
    phi=(cosh(Ek./(2*TTc))).^(-2);
    int=abs(ksi./Ek).^n .* phi *(2*TTc) .* cosh(ksi/(2*TTc)).^2;
end


