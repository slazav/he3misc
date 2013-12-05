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
  pars_cache=pp;

  % read signal
  if 1
    [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile);
  else
    t00 = 3.1;
    f00 = 1512.0;
    df00 = 500.0;
    tf00 = 2000.1;
    a00 = 1.1;
    tx=linspace(0,15, 3000000);
    amp00 = a00*exp(-tx/t00);
    fre00 = f00 + df00*exp(-tx/tf00);
    fre01 = fre00 - df00/tf00*tx.*exp(-tx/tf00);
    xx=amp00.*sin(2*pi*fre00 .* tx);
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
  range_factor=2.1; % integrating range / signal width
  check_factor=3; % large/normal int range
  check_ratio=1.1;  % change of amplitude

  fprintf('running fft (window: %d, step: %d):   0 %%', fft_wind, fft_step);
  j = 1;
  maxi = length(tx)-fft_wind;
  pr1=0;
  bl_wind = blackman(fft_wind)';
  n_level = 0.5;
  time=[]; amp=[]; fre=[]; wid=[];
  df0=[]; fre0=[];
  df1=[]; amp1=[]; fre1=[]; wid1=[]; noise=[];
  pts=[];

  for i=1:fft_step:maxi
    pr2=int32(100*i/maxi);
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); pr1=pr2; end

    [F, A(:,j)] =...
      sigproc.fft(tx(i:i+fft_wind-1), xx(i:i+fft_wind-1).*bl_wind, pp.fmin, pp.fmax);
    time_=tx(round(i+fft_wind/2)); %The recorded times are the middle of each window.
    A(:,j) = abs(A(:,j));

    %
    if pp.fix_lock_in
      A(:,j) = A(:,j)./sigproc.lock_in_gain(F)';
    end


    % make a guess about frequency range (fre0 and df0) 
    % use calculated frequency and width from prev step (fre1 and df1) 
    if (j>2)
      fre0_ = 2*fre1(j-1)-fre1(j-2);
    elseif (j>1)
      fre0_ = fre1(j-1);
    else fre0_ = f0; end

    if ~pp.fixdf && j>1;
      df0_ = df1(j-1);
      if isnan(df0_); break; end
      if pp.decdf && df0_ > df0(j-1); df0_=df0(j-1); end
    else
      df0_ = df;
    end

    if fre0_-df0_ <= pp.fmin ||  fre0_+df0_ >= pp.fmax ;
      fprintf('!!! fmin..fmax range exeeded.\n');
      break;
    end



    % after first pass range is adjust, in second pass actual peak parameters are calculeted

    ii=find(abs(F-fre0_)<df0_);
    [amp1_ fre1_ wid1_ noise_] = find_peak2(F(ii), A(ii,j)');
    df1_=wid1_*range_factor;


    %%% stopping condition:
    if noise_ > amp1_; break; end
    if df1_ >= pp.fmax-pp.fmin; break; end
    if df1_ <= 0; break; end

%    [amp(j), ~, fre(j), wid(j)] = find_peak(F(ii), A(ii,j)');


%    %%% adjast integrating range
%    if (j>1) f0 = 2*fre(j)-fre(j-1);
%    else f0 = fre(j); end
%    if ~pp.fixdf;
%      df1 = wid(j)*range_factor;
%      if isnan(df1); break; end
%      if f0-df1 < pp.fmin df1=f0-pp.fmin; end
%      if f0+df1 > pp.fmax df1=pp.fmax-f0; end
%      if ~pp.decdf || df1<df; df=df1; end
%    end

    %%% find peak (integration in f0 +/- df range)
%    ii=find(abs(F-fre1_)<df1_);
    [amp(j), ampp(j), fre(j), wid(j)] = find_peak(F(ii), A(ii,j)');

%    if ~pp.fixdf;
%      df1 = wid(j)*range_factor;
%      if isnan(df1); break; end
%      if f0-df1 < pp.fmin df1=f0-pp.fmin; end
%      if f0+df1 > pp.fmax df1=pp.fmax-f0; end
%      if ~pp.decdf || df1<df; df=df1; end
%    end

    %%% again adjast integrating range
%    if (j>1) f0 = 2*fre(j)-fre(j-1);
%    else f0 = fre(j); end
%    if ~pp.fixdf; df = wid(j)*range_factor; end
%    if isnan(df); break; end


    fre0(j)=fre0_;
    df0(j)=df0_;
    df1(j)=df1_;
    amp1(j)=amp1_;
    fre1(j)=fre1_;
    wid1(j)=wid1_;
    noise(j)=noise_;
    time(j)=time_;
    pts(j)=length(ii);


%    n_level = peak_noise_level(F, A(:,j), fre0_, df0_*check_factor);
%    n_level = noise_peak_level(F, A(:,j), f0, df*check_factor)
%    if n_level < 1;
%      break;
%    end

    j=j+1;
  end
  fprintf('\b\b\b\b\bok   \n');

  mnoise=mean(noise);
  amp = amp-mnoise*sqrt(pi*pts);

  if plot_steps
    find_figure('tracing steps'); clf; hold on;
%    for j=1:5:length(amp1)
    for j=find(abs(time-40)<0.5)
      ii=find(abs(F-fre0(j))<df0(j));
      plot(F(ii), A(ii,j), '.g-');
      xx = linspace(fre0(j)-df0(j), fre0(j)+df0(j), 1000);
      plot(xx, amp1(j)*exp(-(xx-fre1(j)).^2/wid1(j))+noise(j), 'r-');
    end
  end



  % plot data
  if do_plot
    ff=find_figure('aft'); clf;

    axes(); hold on;
      title([dstr ' ' xfile], 'Interpreter', 'none');
      sigproc.plot_3di(F, time, abs(A)', 'sqrt');
      plot (time, fre1, 'b.-', 'MarkerSize',10);
      plot (time, fre1-wid1/2, 'b-');
      plot (time, fre1+wid1/2, 'b-');
      plot (time, fre1-df0, 'c-');
      plot (time, fre1+df0, 'c-');
      xlim([min(time), max(time)]);
      ylabel('frequency, Hz');
      xlabel('time, s');
      ylim([min(fre)-max(wid), max(fre)+max(wid)]);
    axes('Color', 'none', 'YAxisLocation', 'right'); hold on;

      plot(time, amp1, '.-','LineWidth',2, 'Color', [1,0,0], 'MarkerSize', 16);
      plot(time, amp1, '.-','LineWidth',1, 'Color', [0.2,0,0], 'MarkerSize', 10);

      plot(time, amp, '.-','LineWidth',2, 'Color', [0,1,0], 'MarkerSize', 16);
      plot(time, amp, '.-','LineWidth',1, 'Color', [0,0.2,0], 'MarkerSize', 10);

      xlim([min(time), max(time)]);
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
      figure(ff);
      print('-dpng', '-r72', file_png);
      unix(['chmod 664 ' file_png]);
    end

    save(file_par, 'pars_cache', 'time', 'amp', 'fre', 'wid');
  end

  res.time = time;
  res.fre  = fre;
  res.amp  = amp;
  res.wid  = wid;
  res.pars = pars_cache;
  res.cache = 0;
end




function [amp ampp fre wid] = find_peak(f, a)
  %% Calculate peak moments.
  %% For function amp * exp(-(f-fre).^2/wid^2)
  %% it must find amp, wid, fre.
  %% int is an integral of the peak

  p=abs(a).^2;
%  N = length(f);
%  df = (max(f)-min(f))/N;

  m0 = sum(p);
  m1 = sum(f .* p);
  m2 = sum((f-m1/m0).^2 .* p );

  amp = sqrt(pi*m0);
  fre = m1/m0;
  wid = 2*sqrt(m2/m0);

  ampp = amp/sqrt(wid);
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

function [amp fre wid noise] = find_peak2(f, a)
  % find largest peak (amp, fre), halfwidth (at half height). 
  % error (largest maximum outside peak / peak). 
  aa=abs(a);
  [fmax,amax,imax] = sigproc.max_3pt(f,aa);
  amin=min(aa);
  level=exp(-1);

  % find noise level 
  inr=round((length(aa)-imax)/2);
  in = [1:imax/2 length(aa)-(inr:-1:0)];
  noise = mean(aa(in));

  wid = (sum(aa)/length(f)  ) * (f(end)-f(1)) / sqrt(pi) / (amax-2*noise);

  % recalculate freq
%  ii=find(abs(f-fmax)<wid*3);
%  fre = sum(f(ii).*aa(ii))/sum(aa(ii));

  % smooth signal and recalculate amp and wid 
%  aas=smooth(a(ii), length(ii)/20);
%  amp = max(aas)-2*noise;
%  wid = (sum(aas)/length(f(ii))  ) * (f(ii(end))-f(ii(1))) / sqrt(pi) / amp;

  amp=amax-2*noise;
  fre=fmax;
end
