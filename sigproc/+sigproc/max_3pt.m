function [fm,am,im] = max_3pt(freq, amp)
% find max in quadratic approximation

    [amax,im] = max(amp);

    if im==1 || im==length(freq);
      fm=freq(im); am=amp(im); return;
    end

    % fit amp(freq) near maximum by Ax^2+Bx+C.
    f1 = freq(im-1);
    f2 = freq(im);
    f3 = freq(im+1);
    a1 = amp(im-1);
    a2 = amp(im);
    a3 = amp(im+1);

    AA = ((a2-a1)/(f2-f1)-(a3-a1)/(f3-f1))/...
         ((f2^2-f1^2)/(f2-f1)-(f3^2-f1^2)/(f3-f1));
    BB = (a2-a1)/(f2-f1) - AA*(f2^2-f1^2)/(f2-f1);
    CC = a1 - AA * f1^2 - BB * f1;

    fm = -BB/2.0/AA;
    am = AA*fm^2 + BB*fm + CC;
end
