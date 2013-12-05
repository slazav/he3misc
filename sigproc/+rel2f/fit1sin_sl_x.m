function [time, freq, amp, ph] = fit1sin_sl_x(tx,xx, window, step, f0ran)
% [time(t), amp0(t), freq(2,t), amp(2,t), ph(2,t)] =
%        fit2sin_sl(tx(n),xx(n), window, step, op)
%
% Sliding fit by one sine functions.
% The input is the same as for fft_sl()
%
% -- slazav, feb 2012.

  tx = tx(:);
  xx = xx(:);

  fprintf('running 1sin fit (window: %d, step: %d):   0 %%', window, step);
  maxi = length(tx)-window;
  pr1 = 0;

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

    %[~,p_guess] = fitsin(mean(f0ran),txm,xxm);

    fmfit = fminbnd(@fitsin,f0ran(1),f0ran(2),optimset('Display','off'),...
                    txm,xxm);

    [~,p_fit,par] = fitsin(fmfit,txm,xxm);
    if 0
    find_figure('sin fit');
    clf
    plot(txm,xxm,'ob','MarkerSize',3);
    hold on
    %plot(txm,p_guess,'-k');
    plot(txm,p_fit,'-r');
    pause
    end
    % save values
    time(j) = tx(i+round(window/2));
    amp(j) = hypot(par(1),par(2));
    freq(j) = fmfit;
    
  end
  fprintf('\b\b\b\b\bok   \n');

end

function [resid,val,par] = fitsin(freq,tx,xx)
% first solve for amplitude and phase
  tt = 2*pi*freq*tx;
  A = [sin(tt) cos(tt)];
  par = A \ xx;
  val = A*par;
  resid = sum((val - xx).^2);
end
