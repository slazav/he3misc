function [amp, tau, fit_endt] = rel_fix(time, p2a, t0)
% [amp, tau] = rel_fix(time(n), p2a(n))
%
% fit p2a(:, n) by amp*exp(-t/tau)

    N = length(time);
    tb=time(1);
    te=time(N);
    dt=(te-tb)/(N-1);

    % index for noise averaging (last 2 seconds)
    ii1=round((te-2)/dt)+1;
    if ii1<1 ii1=1; end

    ave_noise = mean(p2a(ii1:N));
    if ave_noise > 0
      fit_endt = time(find(p2a < ave_noise,1));
    else
      fit_endt = time(N);
    end
    fit_ind = find(time >= t0 & time <= fit_endt);
    fit = lsqcurvefit(@fitfun,[1 1],time(fit_ind),p2a(fit_ind));
    tau = fit(2);
    amp = fit(1);
end

%Function for which to make fits of the relaxation times... 
function F = fitfun(x,xdata)
    F = x(1)*exp(-xdata./x(2));
end
