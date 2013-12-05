function F = qfunc(p, x, f0, Q)
    A0=p(1);
    A1=p(2);
    df=f0/Q;
    F = A0 + A1  ./ ( (( x.^2 - f0^2 )/f0/df ).^2 + (x/f0).^2 );% - p(3)*x; % 3rd parameter for linear background
end
