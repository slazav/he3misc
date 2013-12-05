function [p_f, p_a] = find_peaks(freq,amp)
% [p_f(k), p_a(k)] = find_peaks(freq(n),amp(n))
%
% Select local maximums on fft output, find correct frequencies
% and amplitudes using three nearest points quadratic fit.
%
% -- slazav, feb 2012.

  mm=1;
  p_f = 0;
  p_a = 0;

  for (m=2:length(freq)-1)
    % fit amp(freq) near maximum by Ax^2+Bx+C.
    f1 = freq(m-1);
    f2 = freq(m);
    f3 = freq(m+1);

    a1 = amp(m-1);
    a2 = amp(m);
    a3 = amp(m+1);

    if (a2<a1) || (a2<=a3); continue; end

    AA = ((a2-a1)/(f2-f1)-(a3-a1)/(f3-f1))/...
         ((f2^2-f1^2)/(f2-f1)-(f3^2-f1^2)/(f3-f1));
    BB = (a2-a1)/(f2-f1) - AA*(f2^2-f1^2)/(f2-f1);
    CC = a1 - AA * f1^2 - BB * f1;

    f0 = -BB/2.0/AA;
    a0 = AA*f0^2 + BB*f0 + CC;

    p_f(mm) = f0;
    p_a(mm) = a0;
    mm=mm+1;
  end
end
