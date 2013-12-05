function gain = lock_in_gain(fre)
% lock-in filtering function at 100us time constant
% see /rota/Analysis/PS/osc2011/LockIn

  pp = [ -8.6883e-19 3.02828e-16 -4.29819e-14 3.02891e-12 -8.6027e-11 -2.39079e-09...
          3.02678e-07 -1.21324e-05 0.00027199 -0.00374675 0.0331762  -0.248256 1.63534];
  a0 = 5.1312;
  gain = exp(polyval(pp, (fre/1000).^2)) / a0;
  ii=find(abs(fre)>9000);
  gain(ii)=0;
end