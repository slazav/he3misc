function [ps_f, ps_a] = select_side_peaks(p_f, p_a, f0, df,ddf)
% [ps_f(2), ps_a(2)] =
%     select_side_peaks(p_f(n), p_a(n), f0, df, ddf)
%
% Select side peaks at f0+df +/- ddf and f0-df +/- ddf.
%
% -- slazav, feb 2012.

    [~, max1ind] = max(p_a);
    i1=0; a1=-1e99;
    i2=0; a2=-1e99;

    for i=1:length(p_f)
      if (abs(p_f(i)-(f0+df))<ddf)
        if p_a(i) > a1; i1 = i; a1=p_a(i); end
      end
      if (abs(p_f(i)-(f0-df))<ddf)
        if p_a(i) > a2; i2 = i; a2=p_a(i);  end
      end
    end

    if i1 == 0;
      ps_f(1)=f0+df;  ps_a(1)=0;
    else
      ps_f(1)=p_f(i1); ps_a(1)=p_a(i1);
    end

    if i2 == 0;
      ps_f(2)=f0-df;  ps_a(2)=0;
    else
      ps_f(2)=p_f(i2); ps_a(2)=p_a(i2);
    end
end
