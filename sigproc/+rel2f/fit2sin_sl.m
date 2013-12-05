function [time, amp0, freq, amp, ph] = fit2sin_sl(tx,xx, window, step, op)
% [time(t), amp0(t), freq(2,t), amp(2,t), ph(2,t)] =
%        fit2sin_sl(tx(n),xx(n), window, step, op)
%
% Sliding fit by two sine functions.
% The input is the same as for fft_sl()
% Very slow!
% use op value to save/restore results:
%  'save' - save results after fit int 'tmp_sin2fit_sl' file
%  'load' - read previous results instead of fit
%
% -- slazav, feb 2012.

  res_file = 'tmp_sin2fit_sl';

  % save results
  if (op == 'load')
    fprintf('sin2fit: loading previous results\n');
    fo=fopen(res_file, 'r');
    [res] = fscanf(fo, '%g',[8,inf]); 
    for i=1:length(res(1,:))
      time(i) = res(1,i);
      amp0(i) = res(2,i);
      freq(1,i) = res(3,i);
      amp(1,i) = res(4,i);
      ph(1,i) = res(5,i);
      freq(2,i) = res(6,i);
      amp(2,i) = res(7,i);
      ph(2,i) = res(8,i);
    end
    fclose(fo);
    return;
  end


  fprintf('running 2sin fit (window: %d, step: %d):   0 %%', window, step);
  j = 1;
  maxi = length(tx)-window;
  pr1=0;

  % initial values
  [freq0, amp0] = rel2f.fft(tx(1:window), xx(1:window).*blackman(window)');
  [pf,pa] = rel2f.find_peaks(freq0, amp0);
  [p2f, p2a]=rel2f.select_2peaks(pf, pa);

  par(1)=p2a(1)/10; %a0
  par(2)=p2a(1);    %a1
  par(3)=p2a(2);    %a2
  par(4)=1;         %p1
  par(5)=1;         %p2
  par(6)=p2f(1);    %f1
  par(7)=p2f(2);    %f2

  for i=1:step:maxi

    pr2=int32(100*i/maxi);
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); pr1=pr2; end

    txm=tx(i:i+window-1)-tx(i);
    xxm=xx(i:i+window-1);

    par = lsqcurvefit(@func2sin,par,txm,xxm,...
                      [],[],optimset('Display','off'));

    % fix negative amplitude
    if par(2) < 0; par(2) = -par(2); par(4) = par(4) + pi; end
    if par(3) < 0; par(3) = -par(3); par(5) = par(5) + pi; end
    while par(4) > 2*pi; par(4)=par(4)-2*pi; end
    while par(5) > 2*pi; par(5)=par(5)-2*pi; end

    % save values
    time(j) = tx(i);
    amp0(j) = par(1);
    amp(1:2,j) = par(2:3);
    freq(1:2,j) = par(6:7);
    ph(1:2,j) = par(4:5);

    % initial conditions for the next step
    if i+step < length(tx)
      dt=tx(i+step)-tx(i);
      par(4)=par(4) + 2*pi*par(6)*dt; %p1
      par(5)=par(5) + 2*pi*par(7)*dt; %p2
      while par(4) > 2*pi; par(4)=par(4)-2*pi; end
      while par(5) > 2*pi; par(5)=par(5)-2*pi; end
    end
    j=j+1;
  end
  fprintf('\b\b\b\b\bok   \n');

  % save results
  if (op == 'save')
    fprintf('sin2fit: saving results\n');
    fo=fopen(res_file, 'w');
    for i=1:length(time)
      fprintf(fo, '%g  %g %g %g  %g %g %g  %g\n', ...
        time(i), amp0(i), ...
        freq(1,i), amp(1,i), ph(1,i), ...
        freq(2,i), amp(2,i), ph(2,i));
    end
    fclose(fo);
  end

end

% fit function for the fit2sin()
function xx = func2sin(par, tx)
  % par(1) a0
  % par(2) a1
  % par(3) a2
  % par(4) p1
  % par(5) p2
  % par(6) f1
  % par(7) f2
  % par(8) r1
  % par(9) r2
  for i=1:length(tx)
    xx(i) = par(1) + ...
            par(2) *  sin(2*pi*par(6)*tx(i) + par(4)) + ...
            par(3) *  sin(2*pi*par(7)*tx(i) + par(5));
%               x(2) * exp(xdata(i)/x(8)) * sin(2*pi*x(6)*xdata(i) + x(4)) + ...
%               x(3) * exp(xdata(i)/x(9)) * sin(2*pi*x(7)*xdata(i) + x(5));
  end
end
