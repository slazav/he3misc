function [time, freq, amp, ph] = fit1sin_sl(tx,xx, window, step)
% [time(t), amp0(t), freq(2,t), amp(2,t), ph(2,t)] =
%        fit2sin_sl(tx(n),xx(n), window, step, op)
%
% Sliding fit by one sine functions.
% The input is the same as for fft_sl()
%
% -- slazav, feb 2012.

  fprintf('running 1sin fit (window: %d, step: %d):   0 %%', window, step);
  j = 1;
  maxi = length(tx)-window;
  pr1=0;

  % initial values
  [freq0, amp0] = rel2f.fft(tx, xx.*blackman(length(xx))');
  [am, im] = max(amp0);
  fm = freq0(im);

  par(1)=max(abs(xx(1:window)));      %a1
  par(2)=1;         %p1
  par(3)=fm;        %f1
  
  nn = ceil(maxi/step);
  time = zeros(nn,1);
  amp = time;
  freq = time;
  ph = time;
  
  for i=1:step:maxi

    j = (i-1)/step+1;
    pr2=int32(100*i/maxi);
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); pr1=pr2; end

    txm=tx(i:i+window-1)-tx(i);
    xxm=xx(i:i+window-1);

    %p_guess = func1sin(par,txm);

    par = lsqcurvefit(@func1sin,par,txm(:),xxm(:),...
                      [],[],optimset('Display','off','Jacobian','on'));

    if 0
      p_fit = func1sin(par,txm);
      find_figure('sin fit');
      clf
      plot(txm,xxm,'ob','MarkerSize',3);
      hold on
      plot(txm,p_guess,'-k');
      plot(txm,p_fit,'-r');
      pause
    end
    % fix negative amplitude
    if par(1) < 0; par(1) = -par(1); par(2) = par(2) + pi; end
    par(2) = mod(par(2), 2*pi);

    % save values
    time(j) = tx(i+round(window/2));
    amp(j) = par(1);
    freq(j) = par(3);
    ph(j) = mod(par(2) - 2*pi*par(3)*time(j), 2*pi);

    % initial conditions for the next step
    if i+step < length(tx)
      dt=tx(i+step)-tx(i);
      par(2)=par(2) + 2*pi*par(3)*dt;
      par(2) = mod(par(2), 2*pi);
    end
    
  end
  fprintf('\b\b\b\b\bok   \n');

end

% fit function for the fit2sin()
function [F,J] = func1sin(par, tx)
  % par(1) a1
  % par(2) p1
  % par(3) f1
      F1 = sin(2*pi*par(3)*tx + par(2));
      F = par(1) * F1 ;
      if nargout > 1
          F2 = par(1)*cos(2*pi*par(3)*tx + par(2));
          J = [F1 F2 F2*2*pi.*tx];
      end
end
