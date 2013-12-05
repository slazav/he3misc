function res = sig_fit2(dstr, xfile, pars)
% fit relaxation signal

%function [a0, aerr, tau, terr, b, berr, f0, df, tf, ff] = sig_fit(time, fre, amp, varargin)

  % get parameters
  pp.var     = sigproc.par_get('var',     pars, '' );
  pp.t1fit   = sigproc.par_get('t1fit',   pars, 0 );
  pp.t2fit   = sigproc.par_get('t2fit',   pars, pp.t1fit );
  pp.maxa    = sigproc.par_get('maxamp',  pars, 1e9 );
  pp.mina    = sigproc.par_get('minamp',  pars, 0 );
  pp.maxdf   = sigproc.par_get('maxdf',   pars, 0 );
  pp.maxdf_a = sigproc.par_get('maxdf_a', pars, 0 );
  pp.fitfunc = sigproc.par_get('fitfunc', pars, 0 );
  pp.meanfre = sigproc.par_get('meanfre', pars, 0 );
  pp.rembase = sigproc.par_get('rembase', pars, 0 );
  do_plot    = sigproc.par_get('plot',    pars, 1 );
  list_file  = sigproc.par_get('list_file', pars, '' );
  refit      = sigproc.par_get('refit',   pars, 0 );

  % values in pp will go to cache 

  % run sig_trace
  trace_res = sigproc.sig_trace4(dstr, xfile, pars);

  others={};
  if (iscell(xfile))
    others=xfile(2:end);
    xfile=xfile{1};
  end

  % in list_file mode we can use cache!
  if length(list_file)>0
    % set filenames
    pdir = [ list_file '.cache'];
    file_png=[pdir '/fit/' dstr '_' xfile pp.var '.png'];
    file_out=[pdir '/fit/' dstr '_' xfile pp.var '.mat'];

    % check trace result and our cache
    if ~refit && trace_res.cache && unix(['[ -f "' file_png '" -a -s "' file_out '" ]']) == 0
      load(file_out, '-mat', 'pars_cache', 'res');
      if isequal(pp, pars_cache);
%        fprintf('skipping processed file: %s %s\n', dstr, xfile);
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
  int2 = trace_res.int2;
  noise2 = trace_res.noise2;

  if pp.rembase; amp = amp - sqrt(noise2); end;

  if pp.t2fit<=pp.t1fit; pp.t2fit=time(end); end

  % compensate lock-in filter curve -- it must be done in the trace program!
  % amp = amp ./ sigproc.lock_in_gain(fre);

  % find fitting range
  i1=find(time >= pp.t1fit & amp < pp.maxa, 1);
  i2=find(time >= pp.t2fit | amp < pp.mina, 1);
  if length(i2)~=1; i2=length(time); end
  if (pp.maxdf_a>0) i1 = find(fre < min(fre)+pp.maxdf_a,1); end
  ii=i1:i2;  % for amp-time fit

  if (pp.maxdf>0)
    iif1 = find(fre < min(fre)+pp.maxdf);
  else
    iif1 = ii;  % for freq fit
  end
  iif1 = iif1(find(~isnan(fre(iif1))));

  iif=find(time < pp.t2fit & ~isnan(fre));  % for freq-amp fit



%  opts=optimset('Algorithm', 'levenberg-marquardt','DIagnostics', 'on', 'PlotFcns', []);
  if pp.fitfunc==0
%    base = sqrt(mean(noise2(ii)));
    base = amp(end);
    i1 = find((amp(ii)-base)/(amp(ii(1))-base) < 1/exp(1),1);
    if length(i1) && i1>1;  tau = time(ii(i1))-time(ii(1)); else tau=(time(ii(end))-time(ii(1)))/5; end
    par = [max(amp(ii)) * exp(time(ii(1))/tau), tau, base];

    func=@func1exp_b1;
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(func, par, time(ii), amp(ii));
    err = full(sqrt(diag((J'*J)^-1))/length(ii));
    a0  = par(1);
    tau = par(2);
    b   = 0;
    base = abs(par(3));
    aerr = err(1);
    terr = err(2);
    berr = err(3);
  elseif pp.fitfunc==1
    base = amp(end);
    i1 = find((amp(ii)-base)/(amp(ii(1))-base) < 1/exp(1),1);
    if length(i1) && i1>1;  tau = time(ii(i1))-time(ii(1));
    else tau=(time(ii(end))-time(ii(1)))/5; end

    a1=max(amp(ii)) * exp(time(ii(1))/tau);
    par = [a1, tau, base, max(amp)-a1, tau/2];

    ii=1:length(time);
    func=@func1exp_b2;
    [par,resnorm,~,~,~,~,J] =...
       lsqcurvefit(@func1exp_b2, par, time(ii), amp(ii))
    err = full(sqrt(diag((J'*J)^-1))/length(ii));
    a0  = par(1);
    tau = par(2);
    b   = 0;
    base = abs(par(3));
    aerr = err(1);
    terr = err(2);
    berr = err(3);
  end

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
    par3 = [(max(fre(iif1))-min(fre(iif1))) / exp(-2*time(min(iif))/tau), tau/2.0, min(fre(iif1))];

    if (mean(fre(iif1(1:round(end/2)))) < mean(fre(iif1(round(end/2):end))));
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
    xx1=[time(1) time(end)];
    subplot(2,2,1); hold on; title('amp - time');
      plot(time, amp, '*b-', 'MarkerSize',2);
      plot(time(ii), amp(ii), '*r-', 'MarkerSize',2);
      plot(time, func(par,time), 'k-');
      plot(xx1, base*[1 1], 'k--');
      text(0.95, 0.80, txt1, 'Units','normalized', 'HorizontalAlignment','right');
      ylim([0 max(amp)*1.1]);
%      plot(time, amp-func1expt_b(par,time), '*b-', 'MarkerSize',2);
    subplot(2,2,2); hold on; title('freq - time');
      plot(time, f0-fre, '*b-', 'MarkerSize',2);
      plot(time(iif1), f0-fre(iif1), '*m-', 'MarkerSize',2);
      if ~pp.meanfre
        plot(time, f0-func1exp_b(par3,time), 'k-');
        plot(xx1, [0 0], 'k--');
      end
      text(0.95, 0.20, txt3, 'Units','normalized', 'HorizontalAlignment','right');
      ylim([min(f0-fre)-5 max(f0-fre)+5]);
    subplot(2,2,3); hold on; title('log(amp^2), log(fre)');
      plot(time, log(abs(amp.^2-base^2)), '*b-', 'MarkerSize',2);
      plot(time(ii), log(abs(amp(ii).^2-base^2)), '*r-', 'MarkerSize',2);
      plot(time, log(func(par,time).^2-base^2), 'k-');
      plot(time, real(log(fre-f0))-5, '*b-', 'MarkerSize',2);
      plot(time(ii), real(log(fre(ii)-f0))-5, '*m-', 'MarkerSize',2);
      if ~pp.meanfre
        plot(time, real(log(func1exp(par3,time)))-5, 'k-');
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
  res.ampfit=[a0, aerr, tau, terr, base, berr];
  res.frefit=[f0, df, tf];
  res.cache=0;

  if length(list_file)
    if ff>0
      set(ff, 'PaperUnits', 'points', 'PaperSize', [1280,1024],...
                'PaperPosition', [0,0,1280,1024]);
      print('-dpng', '-r72', file_png);
      unix(['chmod 666 ' file_png ' ||:']);
%    set(ff, 'PaperOrientation', 'landscape');
    else
      unix(['touch ' file_png ' ||:']);
      unix(['chmod 666 ' file_png ' ||:']);
    end

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

% exponent
function xx = func1exp1(par, tx)
    xx = exp(-(tx-par(1))/(par(2)));
end

% exponent with base
function xx = func1exp_b(par, tx)
    xx = par(1) *  exp(-tx/(par(2))) + par(3);
end

% exponent with base
function xx = func1exp_b1(par, tx)
    xx = sqrt( (par(1) *  exp(-tx/(par(2)))).^2 + par(3)^2);
end

% exponent with base
function xx = func1exp_b2(par, tx)
    xx = sqrt( (par(1) *  exp(-tx/(par(2))) + par(4) *  exp(-tx/(par(5)))).^2 + par(3)^2);
end

% exponent with base
function xx = func1exp_b1a(par, tx, base)
    xx = sqrt( (par(1) *  exp(-tx/(par(2)))).^2 + base^2);
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
