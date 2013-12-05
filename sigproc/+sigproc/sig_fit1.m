function res = sig_fit1(dstr, xfile, pars)
% fit relaxation signal

%function [a0, aerr, tau, terr, b, berr, f0, df, tf, ff] = sig_fit(time, fre, amp, varargin)

  % get parameters
  pp.t1fit   = sigproc.par_get('t1fit',   pars, 0 );
  pp.t2fit   = sigproc.par_get('t2fit',   pars, pp.t1fit );
  pp.var     = sigproc.par_get('var',     pars, '' );
  pp.maxa    = sigproc.par_get('maxamp',  pars, 1e9 );
  pp.mina    = sigproc.par_get('minamp',  pars, 0 );
  pp.maxdf   = sigproc.par_get('maxdf',  pars, 0 );
  pp.fitfunc = sigproc.par_get('fitfunc', pars, 0 );
  pp.meanfre = sigproc.par_get('meanfre', pars, 0 );
  do_plot    = sigproc.par_get('plot',    pars, 1 );
  list_file  = sigproc.par_get('list_file', pars, '' );
  refit      = sigproc.par_get('refit',   pars, 0 );

  % values in pp will go to cache 


  % run sig_trace
  trace_res = sigproc.sig_trace3(dstr, xfile, pars);

  % in list_file mode we can use cache!
  if length(list_file)>0
    % set filenames
    pdir = [ list_file '.cache'];
    file_png=[pdir '/fit/' dstr '_' xfile pp.var '.png'];
    file_out=[pdir '/fit/' dstr '_' xfile pp.var '.mat'];

    % check trace result and our cache
    if ~refit && trace_res.cache && unix(['[ -s "' file_png '" -a -s "' file_out '" ]']) == 0
      load(file_out, '-mat', 'pars_cache', 'res');
      if isequal(pp, pars_cache);
        fprintf('skipping processed file: %s %s\n', dstr, xfile);
        res.pars=pars; % can be different!
        res.cache=1;
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/fit/']);
    % save original parameters before modifications 
    pars_cache=pp;
  end

  % get data to fit
  time = trace_res.time;
  fre  = trace_res.fre;
  amp  = trace_res.amp;

  if pp.t2fit<=pp.t1fit; pp.t2fit=time(end); end

  % compensate lock-in filter curve -- it must be done in the trace program!
  % amp = amp ./ sigproc.lock_in_gain(fre);

  % find fitting range
  i1=find(time >= pp.t1fit & amp < pp.maxa, 1);
  i2=find(time >= pp.t2fit | amp < pp.mina, 1);
  if length(i2)~=1; i2=length(time); end
  ii=i1:i2;  % for amp-time fit

  if (pp.maxdf>0)
    iif1 = find(fre < min(fre)+pp.maxdf);
  else
    iif1 = ii;  % for freq fit
  end
  iif1 = iif1(find(~isnan(fre(iif1)));

  iif=find(time < pp.t2fit & ~isnan(fre));  % for freq-amp fit

%  opts=optimset('Algorithm', 'levenberg-marquardt','DIagnostics', 'on', 'PlotFcns', []);
  switch pp.fitfunc
  case 0  % exp fit
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1exp, [max(amp),(pp.t2fit-pp.t1fit)/2], time(ii), amp(ii));
    err = full(sqrt(diag((J'*J)^-1))/length(ii));
    err = [err(1), err(2), 0, 0];
    par = [par(1), par(2), 0, 0];
  case 1 %exp fit with base
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1exp_b, [max(amp),(pp.t2fit-pp.t1fit)/5,max(amp)/10], time(ii), amp(ii));
    err = full(sqrt(diag((J'*J)^-1))/length(ii));
    err = [err(1), err(2), 0, err(3)];
    par = [par(1), par(2), 0, par(3)];
  case 2 %exp fit with nonlinear term
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1expt, ...
         [max(amp),(pp.t2fit-pp.t1fit)/5, 0.1], time(ii), amp(ii),...
         [0,0,-inf], [inf,inf,inf]);
    err = full(sqrt(diag((J'*J)^-1))'/length(ii));
    err = [err(1), err(2), err(3), 0];
    par = [par(1), par(2), par(3), 0];
  case 3 %exp fit with nonlinear term and base
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1expt_b, [max(amp),(pp.t2fit-pp.t1fit)/5,0.1,max(amp)/10], time(ii), amp(ii));
    err = full(sqrt(diag((J'*J)^-1))'/length(ii));
  case 4 %fit log(amp with linear func)
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1lin, [1 -1], time(ii), log(amp(ii)));
    err = full(sqrt(diag((J'*J)^-1))'/length(ii));
    err = [err(1), err(2), 0, 0]; %???
    par = [exp(par(1)), -1/par(2), 0, 0];

  case 5 %exp fit with base
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1exp_b, [max(amp)^2, (pp.t2fit-pp.t1fit)/5, max(amp)^2/10], time(ii), amp(ii).^2);
    err = full(sqrt(diag((J'*J)^-1))/length(ii));
    err = [err(1), err(2), 0, err(3)];
    par = [sqrt(par(1)), par(2)/2, 0, sqrt(par(3))];
  end

  a0  = par(1);
  tau = par(2);
  b   = par(3);
  base = par(4);
  aerr = err(1);
  terr = err(2);
  berr = err(3);

  if pp.meanfre
    f00  = mean(fre(iif));
    fa2  = 0;
    fa4  = 0;
    df   = 0;
    tf   = 0;
    f0   = f00;
  else
    par2 = lsqcurvefit(@func1a2f, ...
      [min(fre(iif)), a0/(max(fre(iif))-min(fre(iif))), 0.1], ...
      fre(iif), amp(iif).^2);

    f00  = par2(1);
    fa2  = par2(2);
    fa4  = par2(3);
    %  par3 = lsqcurvefit(@func1expf, [min(fre),max(fre)-min(fre),tau/2.0], time(iif1), fre(iif1));
    % fit freq with exp with fixed f0!


    par3 = [(max(fre(iif1))-min(fre(iif1))) / exp(-2*time(min(iif))/tau), tau/2.0, min(fre(iif1))]

    if (mean(fre(iif1(1:end/2))) < mean(fre(iif1(end/2:end))));
       par3(1)=-par3(1); par3(3)=max(fre(iif1));
    end

    par3 = lsqcurvefit(@func1exp_b, par3, time(iif1), fre(iif1));
    df   = par3(1);
    tf   = par3(2);
    f0   = par3(3);
  end


%  tau_window=5;

%  for t=time(1),tau_step,time(end)
%    ii=find(time > t-tau_window/2 & time <= pp.t2fit);
%    [par,resnorm,~,~,~,~,J] =...
%       lsqcurvefit(@func1expt, [100,15,0.2,0.1], time(ii), amp(ii));

%    p = lsqcurvefit
%    (@func1expf, [min(fre),max(fre)-min(fre),tau/2.0], time(iif), fre(iif));
%  end

  if do_plot==1

    txt1={
      'Amp'
      ['Amp0:     ' num2str(a0) ' +/- ' num2str(aerr) ]
      ['Tau:      ' num2str(tau) ' +/- ' num2str(terr) ]
      ['nonexp:   ' num2str(b) ' +/- ' num2str(berr) ]
      ['Base:     ' num2str(base) ]
    };
    txt2={
      'Fre'
      ['F0:     ' num2str(f00) ]
      ['dF/dA2: ' num2str(fa2) ]
      ['dFdA4:  ' num2str(fa4) ]
    };
    txt3={
      'Fre'
      ['F0:     ' num2str(f0) ]
      ['dF:     ' num2str(df) ]
      ['tau_f:  ' num2str(tf) ]
    };

    ff=find_figure('af'); clf; hold on;
    subplot(2,2,1); hold on; title('amp - time');
      plot(time, amp, '*b-', 'MarkerSize',2);
      plot(time(ii), amp(ii), '*r-', 'MarkerSize',2);
      plot(time, func1expt_b(par,time), 'k-');
      plot(time, func1exp(par,time), 'k--');
%      plot(time, amp-func1expt_b(par,time), '*b-', 'MarkerSize',2);
      text(0.95, 0.80, txt1, 'Units','normalized', 'HorizontalAlignment','right');
      ylim([0 max(amp)*1.1]);
    subplot(2,2,2); hold on; title('freq - time');
      plot(time, f0-fre, '*b-', 'MarkerSize',2);
      plot(time(iif1), f0-fre(iif1), '*m-', 'MarkerSize',2);
      if ~pp.meanfre
        plot(time, f0-func1exp_b(par3,time), 'k-');
      end
      text(0.95, 0.20, txt3, 'Units','normalized', 'HorizontalAlignment','right');
      ylim([min(f0-fre)-5 max(f0-fre)+5]);
    subplot(2,2,3); hold on; title('log(amp) - time');
      plot(time, log(amp-base), '*b-', 'MarkerSize',2);
      plot(time(ii), log(amp(ii)-base), '*r-', 'MarkerSize',2);
      plot(time, log(func1expt_b(par,time)-base), 'k-');
      plot(time, log(func1exp(par,time)-base), 'k--');
      plot(time, real(log(fre-f0))/2-5, '*b-', 'MarkerSize',2);
      plot(time(ii), real(log(fre(ii)-f0))/2-5, '*m-', 'MarkerSize',2);
      if ~pp.meanfre
        plot(time, real(log(func1exp(par3,time)))/2-5, 'k-');
      end
      text(0, 1.2, [dstr ' ' xfile], 'Units', 'normalized', 'interpreter', 'none');
%      ylim([log(min(amp-base))-0.2 log(max(amp-base))+0.5]);
    subplot(2,2,4); hold on; title('amp^2(freq)');
      xx=linspace(min(fre), max(fre)+10, 100);
      plot(fre, amp.^2, '*b-', 'MarkerSize',2);
      if ~pp.meanfre
        plot(xx, func1a2f(par2,xx), 'k-');
      end
      text(0.05, 0.80, txt2, 'Units','normalized');
      xlim([min(xx) max(xx)]);
      ylim([0 max(amp)^2 * 1.1]);
  else
    ff=0
  end

  res.dstr=dstr;
  res.file=xfile;
  res.pars=pars;
  res.ampfit=[a0, aerr, tau, terr, b, berr];
  res.frefit=[f0, df, tf];
  res.cache=0;

  if length(list_file)
    set(ff, 'PaperUnits', 'points', 'PaperSize', [1280,1024],...
              'PaperPosition', [0,0,1280,1024]);
    print('-dpng', '-r72', file_png);
    unix(['chmod 666 ' file_png ' ||:']);
%    set(ff, 'PaperOrientation', 'landscape');

    save(file_out, 'pars_cache', 'res');
  end


end


% exponent for freq fit. 
function xx = func1expf(par, tx)
%    il = find(tx<=0);
%    ir = find(tx>0);
%    xx(il) = par(1) + par(2);
%    xx(ir) = par(1) + par(2) *  exp(-tx(ir)/(par(3)));
    xx = par(1) + par(2) *  exp(-tx/(par(3)));
end

% function for fit f-a dep.
function f = func1fa2(par, a2)
    f = par(1) + par(2) * a2 + par(3)*a2.^2;
end

function a2 = func1a2f(par, f)
    f0=par(1);
    a2 = par(2) * (f-f0) + par(3)*(f-f0).^2;
end


% exponent 
function xx = func1exp(par, tx)
    xx = par(1) *  exp(-tx/(par(2)));
end

% exponent with base
function xx = func1exp_b(par, tx)
    xx = par(1) *  exp(-tx/(par(2))) + par(3);
end

% exponent with relaxation time dependance 
function xx = func1expt(par, tx)
  A=par(1);
  tau=par(2);
  b=par(3);
  xx=zeros(size(tx));
  ii=find(exp(2*tx/tau) > b*A^2);
  xx(ii)=A ./ sqrt( exp(2*tx(ii)/tau) - b*A^2);
end

% exponent with relaxation time dependance and base
function xx = func1expt_b(par, tx)
  A=par(1);
  tau=par(2);
  b=par(3);
  A0=par(4);
  xx=zeros(size(tx));
  ii=find(exp(2*tx/tau) > b*A^2);
  xx(ii) = A ./ sqrt(exp(2*tx(ii)/tau) - b*A^2) + A0;
end

function xx = func1lin(par, tx)
    xx = par(1) + par(2)*tx;
end
