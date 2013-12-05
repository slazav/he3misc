function [fre, amp] = fft_prec(t, x, f1, f2, df)
% fft_prec(t, x, f1, f2, df)
%
% Do fourier transformation of x(t) in range f1..f2 with step df.

    j=1;
    for f = f1 : df : f2
      amp(j) = sum(x .* (cos(2.0*pi*f*t) - sin(2.0*pi*f*t) * i)) * 2.0/(length(t)-1);
      fre(j) = f;
      j=j+1;
    end
end


