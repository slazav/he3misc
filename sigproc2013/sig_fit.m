function sig_fit(varargin)
%% Trace relaxation signals.
  func=@sig_fit03;
  sigproc2013.sig_process_list(func, 'readonly=0', varargin{:});
end


function data = sig_fit03(data, list_file)

  parstr=data.pars;

  pars = sigproc2013.sig_fit03_pars(parstr); % read+fft3d+trace parameters
  % our own parameters

  % other parameters are not saved in the cache
  refit          = sigproc2013.par_get('refit',       data.pars, 0 ); % force processing
  retrace        = sigproc2013.par_get('retrace',     data.pars, 0 ); % force processing
  do_plot        = sigproc2013.par_get('do_plot',     data.pars, 1 ); % plot picture
  do_png         = sigproc2013.par_get('do_png',      data.pars, 1 ); % create png files

  % wi need refit if retrace is done
  if retrace; refit=1; end

  % run sig_trace for one entry
  data = sigproc2013.sig_trace05(data, list_file);

  % refit if there are no fit field
  % or parameters are different

  if ~isfield(data, 'fit') ||...
     ~isfield(data.fit, 'pars') || ~isequal(data.fit.pars, pars)
     refit=1;
  end

  % do refit if we need png picture
  if length(list_file) && do_plot && do_png
    pdir = [ list_file '.cache/fit'];
    unix(['mkdir -p -m 775 -- ' pdir ]);
    file_png=[pdir '/' data.alias '.png'];
    if unix(sprintf('[ -s "%s" ]', file_png)) ~= 0; refit=1; end
  end

  if ~refit; return; end;

  % save original parameters before modifications
  res.pars  = pars;

  % get data to fit
  time   = data.trace.time;
  freq   = data.trace.freq;
  amp2   = data.trace.amp2;
  noise2 = data.trace.noise2;

  % set fitting ranges
  if pars.t2fit <= pars.t1fit; pars.t2fit = time(end); end
  if pars.t2fre <= pars.t1fre; pars.t2fre = time(end); end
  if pars.t2af  <= pars.t1af;  pars.t2af  = time(end); end

  i1a=find(~isnan(amp2) & time >= pars.t1fit & sqrt(amp2) < pars.maxamp, 1);
  i2a=find( isnan(amp2) | time >= pars.t2fit | sqrt(amp2) < pars.minamp, 1);
  i1f=find(~isnan(freq) & time >= pars.t1fre & sqrt(amp2) < pars.maxamp_fre, 1);
  i2f=find( isnan(freq) | time >= pars.t2fre | sqrt(amp2) < pars.minamp_fre, 1);
  i1x=find(~isnan(freq) & ~isnan(amp2) & time >= pars.t1af & sqrt(amp2) < pars.maxamp_af, 1);
  i2x=find( isnan(freq) |  isnan(amp2) | time >= pars.t2af | sqrt(amp2) < pars.minamp_af, 1);

  if length(i1a)~=1; i1a=1; end
  if length(i1f)~=1; i1f=1; end
  if length(i1x)~=1; i1x=1; end
  if length(i2a)~=1; i2a=length(time); end
  if length(i2f)~=1; i2f=length(time); end
  if length(i2x)~=1; i2x=length(time); end

  fit_opts=statset('Display', 'off');

  %%% frequency fit
  iif = i1f:i2f-1;
  if pars.maxdf_fre > 0; iif = iif(find(freq(iif) <= min(freq(iif))+pars.maxdf_fre)); end
  if pars.maxdf_fre < 0; iif = iif(find(freq(iif) >= max(freq(iif))+pars.maxdf_fre)); end
  if length(iif)<2;  error('Empty fit range for frequency.'); end
  if pars.func_fre==0 % mean freq
    res.f0     = mean(freq(iif));
    res.f0_err = sqrt( sum( (freq(iif)-res.f0).^2 )/length(iif) );
    res.df     = 0;
    res.df_err = 0;
    res.tf     = 0;
    res.tf_err = 0;
    % make standard function for plotting
    ppf=[res.f0];
    func_fre = @(p,x)( p(1)*ones(size(x)));
  elseif pars.func_fre==1 % expotent fit
    % function
    func_fre = @(p,x)( p(1) + p(2)*exp(-x/p(3)) );
    % initial conditions
    t=(time(iif(end))-time(iif(1)))/3;
    df=freq(iif(1))-freq(iif(end));
    ppf=[freq(iif(end))-14, df*exp(time(iif(1))/t), t];
    % fit
    [ppf,rr,J,~,~] = nlinfit(time(iif),freq(iif),func_fre,ppf,fit_opts);
    ci = nlparci(ppf,rr,'Jacobian',J);
    errf = (ci(:,2)-ci(:,1))'/2; %'
    res.f0     = ppf(1);
    res.f0_err = errf(1);
    res.df     = ppf(2);
    res.df_err = errf(2);
    res.tf     = ppf(3);
    res.tf_err = errf(3);
  else
    error('Unknown func_fre value.');
  end

  %%%% amplitude fit
  iia = i1a:i2a;
  if pars.maxdf_amp > 0; iia = iia(find(freq(iia) <= res.f0+pars.maxdf_amp,1):end); end
  if pars.maxdf_amp < 0; iia = iia(find(freq(iia) >= res.f0+pars.maxdf_amp,1):end); end
  if length(iia)<2;  error('Empty fit range for amplitude.'); end

  anoise2=mean(noise2(iia));
  if pars.fixnoise; amp2a = amp2(iia)-noise2(iia);
  else              amp2a = amp2(iia)-anoise2;
  end

  if pars.func_amp==0 % fit amp2-mean(noise2) or amp2-noise2 with exp without base
    % function
    func_amp = @(p,x)( p(1)*exp(-x/p(2)) );
    % initial conditions
    t=(time(iia(end))-time(iia(1)))/6;
    da=amp2(iia(1))-amp2(iia(end));
    ppa=[da*exp(time(iia(1))/t), t];
    [ppa,rr,J,~,~] = nlinfit(time(iia), amp2a,func_amp,ppa, fit_opts);
    ci = nlparci(ppa,rr,'Jacobian',J);
    erra = (ci(:,2)-ci(:,1))'/2; %'
    res.amp      = sqrt(ppa(1));
    res.amp_err  = erra(1);
    res.tau      = ppa(2)*2;
    res.tau_err  = erra(2)*2;
    res.base     = sqrt(anoise2);
    res.base_err = sqrt( sum( (noise2(iia)-res.base).^2 )/length(iia) );

  elseif pars.func_amp==2 % 2-exp fit of amp2-mean(noise2) or amp2-noise2 without base
    % function
    func_amp = @(p,x)( p(1)*exp(-x/p(2)) + p(3)*exp(-x/p(4)));
    % initial conditions
    if pars.tcrit >= time(iia(end)) ||...
       pars.tcrit <= time(iia(1)); pars.tcrit=mean(time(iia([1,end]))); end
    icrit=find(time>pars.tcrit,1);
    if length(icrit)==0; icrit=iia(round(end/2)); end

    t1=(time(iia(end))-time(icrit))/3;
    t2=(time(icrit)-time(iia(1)))/3;
    da1=(amp2(icrit)-amp2(iia(end)) ) * exp(time(icrit)/t1);
    da2=(amp2(iia(1))-amp2(icrit)) * exp(time(iia(1))/t2) - da1;
    ppa=[da1 t1 da2 t2];
    erra = [0 0 0 0];
    if 1
      [ppa,rr,J,~,~] = nlinfit(time(iia), amp2a,func_amp,ppa, fit_opts);
      ci = nlparci(ppa,rr,'Jacobian',J);
      erra = (ci(:,2)-ci(:,1))'/2; %'
    end
    res.amp      = sqrt(ppa(1));
    res.amp_err  = erra(1);
    res.tau      = ppa(2)*2;
    res.tau_err  = erra(2)*2;
    res.base     = sqrt(anoise2);
    res.base_err = sqrt( sum( (noise2(iia)-res.base).^2 )/length(iia) );
    res.amp2      = sqrt(ppa(3));
    res.amp2_err  = erra(3);
    res.tau2      = ppa(4)*2;
    res.tau2_err  = erra(4)*2;
    res.t_crit  = (res.tau*res.tau2)/(res.tau-res.tau2) * log(res.amp2/res.amp * 4);
    res.a2_crit = func_amp(ppa, res.t_crit);
    res.df_crit  = res.f0-interp1(time, freq, res.t_crit);
  else
    error('Unknown func_amp value.');
  end

  %%%% amp-fre fit
  iix = i1x:i2x;
  if pars.maxdf_af > 0; iix = iix(find(freq(iix) <= res.f0+pars.maxdf_af,1):end); end
  if pars.maxdf_af < 0; iix = iix(find(freq(iix) >= res.f0+pars.maxdf_af,1):end); end
  if length(iix)<2;  error('Empty fit range for freq vs amplitude.'); end
  if pars.func_af==0 %
    % function
    func_af = @(p,x)( p(2)*(x-p(1)) + p(3)*(x-p(1)).^2 );
    ppx=[res.f0, 1e-5, 1e-5];
    [ppx,rr,J,~,~] = nlinfit(freq(iix), amp2(iix),func_af, ppx, fit_opts);
    ci = nlparci(ppx,rr,'Jacobian',J);
    errx = (ci(:,2)-ci(:,1))'/2; %'
    res.af_f0      = ppx(1);
    res.af_f0_err  = errx(1);
    res.af_a1      = ppx(2);
    res.af_a1_err  = errx(2);
    res.af_a2      = ppx(3);
    res.af_a2_err  = errx(3);
  else
    error('Unknown func_af value.');
  end


  %%%% plot a picture

  if do_plot==1
    txt_amp={
      sprintf('A0:  %.4f +/- %.4f', res.amp, res.amp_err)
      sprintf('tau: %.4f +/- %.4f s', res.tau, res.tau_err)
      sprintf('base: %.4f +/- %.4f', res.base, res.base_err)
    };
    txt_af={
      sprintf('F0:  %.4f +/- %.4f Hz', res.af_f0, res.af_f0_err)
      sprintf('A1:  %.2e +/- %.2e', res.af_a1, res.af_a1_err)
      sprintf('A2:  %.2e +/- %.2e', res.af_a2, res.af_a2_err)
    };
    txt_fre={
      sprintf('F0: %.2f +/- %.2f Hz', res.f0, res.f0_err)
      sprintf('dF: %.2f +/- %.2f Hz', res.df, res.df_err)
      sprintf('tau_f: %.2f +/- %.2f s ', res.tf, res.tf_err)
    };
    ff=find_figure('af'); clf; hold on;
    tt=[time(1) time(end)];
    subplot(3,2,1); hold on; ylabel('amp vs time');
      plot(time,      sqrt(amp2), '*b-', 'MarkerSize',2);
      plot(time(iia), sqrt(amp2(iia)), '*r-', 'MarkerSize',2);
      plot(time, sqrt(func_amp(ppa,time) + anoise2), 'k-');
      plot(tt, res.base*[1 1], 'k--');
      text(0.95, 1.0, txt_amp, 'Units','normalized',...
        'HorizontalAlignment','right', 'VerticalAlignment','top');
      ylim([0 sqrt(max(amp2))*1.1]);
      xlim(tt)
      if pars.func_amp==2;
        plot(time, sqrt(ppa(1)*exp(-time/ppa(2)) + anoise2), 'k--');
        plot(res.t_crit, sqrt(res.a2_crit + anoise2), '*g');
      end
    subplot(3,2,2); hold on; ylabel('freq vs time');
      plot(time,      res.f0-freq, '*b-', 'MarkerSize',2);
      plot(time(iif), res.f0-freq(iif), '*m-', 'MarkerSize',2);
      plot(time,      res.f0-func_fre(ppf, time), 'k-');
      text(0.95, 0.00, txt_fre, 'Units','normalized',...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom');
      ylim([min(res.f0-freq)-5 max(res.f0-freq)+5]);
      xlim(tt)
      if pars.func_amp==2;
        plot(res.t_crit, res.df_crit, '*g');
      end

    subplot(3,2,3); hold on; ylabel('\delta amp2 vs time');
      plot(time,      amp2 - func_amp(ppa,time) , '*b-', 'MarkerSize',2);
      plot(time(iia), amp2(iia) - func_amp(ppa,time(iia)), '*r-', 'MarkerSize',2);
      plot(min(time), amp2(iia) - func_amp(ppa,time(iia)), '*r-', 'MarkerSize',2);
      plot(tt,[0 0], 'k-');
      xlim(tt)

    subplot(3,2,4); hold on; ylabel('\delta freq,Hz vs time');
      plot(time,      freq - func_fre(ppf, time), '*b-', 'MarkerSize',2);
      plot(time(iif), freq(iif) - func_fre(ppf, time(iif)), '*m-', 'MarkerSize',2);
      plot(tt,[0 0], 'k-');
      xlim(tt)

    subplot(3,2,5); hold on; ylabel('log(amp^2), log(fre)');
      plot(time,      log(amp2), '*b-', 'MarkerSize',2);
      plot(time(iia), log(amp2(iia)), '*r-', 'MarkerSize',2);
      plot(time,      log(func_amp(ppa,time)+anoise2), 'k-');
      plot(time,      real(log(freq-res.f0))-5, '*b-', 'MarkerSize',2);
      plot(time(iif), real(log(freq(iif)-res.f0))-5, '*m-', 'MarkerSize',2);
      plot(time,      real(log(func_fre(ppf,time)-res.f0))-5, 'k-');
      ylim([log(min(amp2)) log(max([amp2 res.df * exp(-5)]))]);
      xlim(tt)
      if pars.func_amp==2;
        plot(res.t_crit, log(res.a2_crit), '*g');
      end
    subplot(3,2,6); hold on; ylabel('amp^2(freq)');
      xx=linspace(min(freq), max(freq)+10, 100);
      plot(freq, amp2, '*b-', 'MarkerSize',2);
      plot(freq(iix), amp2(iix), '*r-', 'MarkerSize',2);
      plot(freq, func_af(ppx,freq), 'k-');
      text(0.05, 0.80, txt_af, 'Units','normalized');
      xlim([min(freq) max(freq)])
      ylim([0 max(amp2) * 1.1]);

    subplot(3,2,2);
      if iscell(data.file); f=data.file{1}; else f=data.file; end
      f={data.alias; f};
      text(0,1, f, 'units','normalized',...
         'Interpreter','none', 'HorizontalAlignment', 'center',...
         'VerticalAlignment', 'bottom');
  else
    ff=0
  end

  %%%% create png if needed
  if length(list_file) && ff>0 && do_png
    set(ff, 'PaperUnits', 'points', 'PaperSize', [800,600],...
            'PaperPosition', [0,0,800,600]);
    print('-dpng', '-r72', file_png);
    unix(['chmod 664 ' file_png]);
  end

  data.fit = res;
end

