function res = sig_trace1(dstr, xfile, pars)
% Trace sliding fft signals.
% Parameters:
%  f0     - starting frequency (default: auto, right of two highest peaks)
%  t1 t2  - time range (default: 0..end of signal)
%  df     - tracing freq window (default: 50Hz)
%  dfs    - tracing freq window (default: 10Hz)
%  window - sliding fft window (default: 1/100 of signal length)
%           (fft step is 1/10 of window)
%  remconst - remove constant frequency peaks (0 or 1, default: 0)
%  plot     - plot figure with results (0 or 1, default: 1)
%  list_file - put result in the <list_file>.trace/<dstr>_<xfile>.txt .png
%
%  res.cache is 1 or 0 depending on if cache was used or not

  % in the list mode we can use cache
  % do nothing if txt, png, par files existst and pars are same
  % get parameters:

  pp.f0    = sigproc.par_get('f0',  pars, 0 );
  pp.fmin  = sigproc.par_get('fmin',  pars, 10 );
  pp.fmax  = sigproc.par_get('fmax',  pars, 3000 );
  pp.df  = sigproc.par_get('df',  pars, 0 );
  pp.t0  = sigproc.par_get('t0',  pars, 0 );
  pp.t1  = sigproc.par_get('t1',  pars, 0 );
  pp.t2  = sigproc.par_get('t2',  pars, pp.t1 );
  pp.fix_lock_in = sigproc.par_get('fix_lock_in', pars, 1 ); % fix amplitude according to lock-in curve
  pp.fixdf       = sigproc.par_get('fixdf',       pars, 0 ); % use constant df as initial guess for signal width
  pp.decdf       = sigproc.par_get('decdf',       pars, 0 ); % do not increase width guess
  pp.var         = sigproc.par_get('var',         pars, ''); % variant name (if there are several signals in one)
  pp.window      = sigproc.par_get('window',      pars, 0 );
  % other parameterss are not saved in the cache
  do_plot        = sigproc.par_get('plot',        pars, 1 );
  list_file      = sigproc.par_get('list_file',   pars, '');
  refit          = sigproc.par_get('refit',       pars, 0 ); % don't use cache
  usecache       = sigproc.par_get('usecache',    pars, 0 ); % always use cache
  plot_steps     = sigproc.par_get('plot_steps',  pars, 0 ); %
  model          = sigproc.par_get('model',       pars, 0 ); %

  % get data from cache if possible
  if length(list_file)
    pdir = [ list_file '.cache'];
    file_png=[pdir '/trace/' dstr '_' xfile pp.var '.png'];
    file_par=[pdir '/trace/' dstr '_' xfile pp.var '.mat'];
    if ~refit && unix(['[ -s "' file_png '" -a -s "' file_par '" ]']) == 0
      load(file_par, '-mat', 'pars_cache', 'time', 'fre', 'amp', 'wid');
      if usecache || isequal(pp, pars_cache)
        fprintf('skipping processed file: %s %s\n', dstr, xfile);
        res.time = time;
        res.fre  = fre;
        res.amp  = amp;
        res.wid  = wid;
        res.pars = pars_cache;
        res.cache = 1;
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/trace/']);
  end
  % save original parameters before modifications
  time=[]; amp=[]; fre=[]; wid=[];
  pars_cache=pp;

  % read signal
  if model==0
    [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile, pars);
  else
    t00 = 3.1;
    f00 = 1512.0;
    df00 = 500.0;
    tf00 = 1000;
    a00 = 1.1;
    a00n = 3;
    tx=linspace(0,15, 3000000);
    dt=(tx(end)-tx(1))/length(tx);
    amp00 = a00*exp(-tx/t00);
    fre00 = f00 + df00*exp(-tx/tf00);

    ph00 = 2*pi*( tx*f00 - df00*tf00*exp(-tx/tf00) );
    noise = random('normal', 0, a00n, size(tx));
    xx=amp00.*sin(ph00) + noise;

%oo=fopen('test.dat', 'w');
%fprintf(oo, '%e %e %e %e %e\n', [tx; xx; fre; amp; noise]);
%fprintf(oo, '%d\n', dt);
%fprintf(oo, '%d\n', round(xx*1000));
%fclose(oo);
 end


  % set default t1, t2 and window
  if pp.t2<=pp.t1; pp.t2=tx(end); end
  if pp.window<=0; pp.window=floor(length(tx)/100); end

  % subtract t0, cut to [t1 .. t2]
  tx=tx-pp.t0;
  ii=find(tx>=pp.t1 & tx<=pp.t2);
  tx=tx(ii); xx=xx(ii);

  % find starting frequency and df if it is not set yet
  if pp.f0==0 || pp.df==0
    [f, a] = sigproc.fft(tx(1:pp.window), xx(1:pp.window),pp.fmin,pp.fmax);
%    sa=smooth(abs(a),length(a)/50);
    sa=a;
    [~, f0, df] = find_peak1(f, sa, 0.2);
    if pp.df==0; pp.df=df*2; end
    if pp.f0==0; pp.f0=f0; end
  end

  % print parameters
  fprintf('tracing parameters: %s\n', pars);

  % do sliding fft + tracing
  fft_wind = pp.window;
  fft_step=fft_wind/10;
  f0=pp.f0;
  df=pp.df;

  time=[]; amp=[]; fre=[]; noise=[];

  fprintf('running fft (window: %d, step: %d):   0 %%', fft_wind, fft_step);
  j = 1;
  maxi = length(tx)-fft_wind;
  pr1=0;
  bl_wind = blackman(fft_wind)';
  n_level = 0.5;

  for i=1:fft_step:maxi
    pr2=int32(100*i/maxi);
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); pr1=pr2; end

    % do fft
    [F, A(:,j)] =...
      sigproc.fft(tx(i:i+fft_wind-1), xx(i:i+fft_wind-1).*bl_wind, pp.fmin, pp.fmax);
    time(j)=tx(round(i+fft_wind/2)); %The recorded times are the middle of each window.
    A(:,j) = abs(A(:,j));

    % apply lock-in correction
    if pp.fix_lock_in; A(:,j) = A(:,j)./sigproc.lock_in_gain(F)'; end

    % after first pass range is adjust, in second pass actual peak parameters are calculeted
    ii=find(abs(F-f0)<df);
    [~, fre(j), noise(j), snr] = find_peak2(F(ii), A(ii,j)');
    f0=fre(j);

    %%% stopping condition:
%    if snr < 2; break; end
%    if j>1 && abs(fre(j)-fre(j-1))>df/5; break; end

    j=j+1;
  end
  fprintf('\b\b\b\b\bok   \n');


  % second pass -- integrate signal
  mnoise=mean(noise); %mean noise
  atime=time;
  for j=1:length(fre)
    %%% find peak (integration in f0 +/- df range)
    ii=find(abs(F-fre(j))<df);
    amp(j)   = sqrt(sum(abs(A(ii,j)).^2));
    noise(j) = sqrt(length(ii)*mnoise);
    if (amp(j)/noise(j)<1)
      amp=amp(1:j-1);
      fre=fre(1:j-1);
      atime=time(1:j-1);
      noise=noise(1:j-1);
      break;
    end
%    amp(j) = amp(j) - noise(j);
  end

  if plot_steps
    find_figure('tracing steps'); clf; hold on;
    plot(F, A(:,j), '.r-');
    mna=min(A(:,j));
    mxa=max(A(:,j));
    plot(f0+[-1 -1 1 1]*df, [mna mxa mxa mna], '.b-');
  end


  % plot data
  if do_plot
    ff=find_figure('aft'); clf;

    axes(); hold on;
      title([dstr ' ' xfile], 'Interpreter', 'none');
      sigproc.plot_3di(F, time, abs(A)', 'sqrt');
      plot (atime, fre, 'b.-', 'MarkerSize',10);
      plot (atime, fre-df, 'b-');
      plot (atime, fre+df, 'b-');
      xlim([min(time), max(time)]);
      ylabel('frequency, Hz');
      xlabel('time, s');
      ylim([min(fre)-df, max(fre)+df]);
    axes('Color', 'none', 'YAxisLocation', 'right'); hold on;
      plot(atime, amp, '.-','LineWidth',2, 'Color', [0,1,0], 'MarkerSize', 16);
      plot(atime, amp, '.-','LineWidth',1, 'Color', [0,0.2,0], 'MarkerSize', 10);
      plot(atime, noise, '.-','LineWidth',2, 'Color', [1,1,0], 'MarkerSize', 16);
      plot(atime, noise, '.-','LineWidth',1, 'Color', [0.2,0.2,0], 'MarkerSize', 10);
      xlim([min(time), max(time)]);
      ylim([0, max(amp)]);
      ylabel('amplitude, V');
      xlabel('time, s');
  else ff=0
  end

  % put txt and png into cache in the list mode
  if length(list_file)
    % save text file with time-fre-amp.
%    fo=fopen(file_txt, 'w');
%    fprintf(fo, '# %s %s\n', dstr, xfile);.
%    fprintf(fo, '# done by sig_trace.m program\n');.
%    fprintf(fo, '# parameters: %s\n', pars);.
%    fprintf(fo, '#%10s %10s %10s\n', 'time', 'freq', 'amp');.
%    for j=1:length(time)
%      fprintf(fo, '%f %f %f\n', time(j), fre(j), amp(j));.
%    end
%    fclose(fo);

    % save png file.
    if ff>0
      set(ff, 'PaperUnits', 'points', 'PaperSize', [1280,1024],...
              'PaperPosition', [0,0,1280,1024]);
      set(0,'CurrentFigure',ff);
      print('-dpng', '-r72', file_png);
      unix(['chmod 664 ' file_png]);
    end

    save(file_par, 'pars_cache', 'time', 'amp', 'fre', 'wid');
  end

  res.time = atime;
  res.fre  = fre;
  res.amp  = amp;
  res.wid  = wid;
  res.pars = pars_cache;
  res.cache = 0;
end

function [amp fre wid] = find_peak1(f, a, level)
  % find largest peak (amp, fre), halfwidth (at half height)
  % error (largest maximum outside peak / peak)
  aa=abs(a);
  [fmax,amax,imax] = sigproc.max_3pt(f,aa);
  amin=min(aa);
  % find points closiest to the maximum where (a-amin) < level (amax-amin)
  i1 = find(aa(imax:-1:1)-amin < (amax-amin)*level, 1);
  if length(i1); i1=imax-i1+1; else i1=1; end
  i2 = find(aa(imax:end)-amin  < (amax-amin)*level, 1);
  if length(i2); i2=imax+i2-1; else i2=length(a); end
  amp=amax-amin;
  fre=fmax;
  wid=(f(i2)-f(i1))/2;
end

function [amp fre noise snr] = find_peak2(f, a)
  [amp,imax] = max(a);
  fre=f(imax);

  ir = round((length(a)-imax)/2);
  il = round(imax/2);
  iin = [1:il, length(a)-(ir:0)]; % noise
  iis = [il:length(a)-ir];        % signal
  noise  = mean(a(iin).^2);
  signal = mean(a(iis).^2);
  snr   = signal/noise;
end


