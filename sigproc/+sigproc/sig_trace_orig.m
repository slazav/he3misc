function res = sig_trace1(dstr, xfile, pars)
% Trace sliding fft signals.
% "Original" version, where output is in matrix format.
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

  % in the list mode we can use cache
  % do nothing if txt, png, par files existst and pars are same
  % get parameters:
  pp.f0  = sigproc.par_get('f0',  pars, 0 );
  pp.fmin  = sigproc.par_get('fmin',  pars, 10 );
  pp.fmax  = sigproc.par_get('fmax',  pars, 3000 );
  pp.df  = sigproc.par_get('df',  pars, 50 );
  pp.dfs = sigproc.par_get('dfs', pars, 10 );
  pp.t0  = sigproc.par_get('t0',  pars, 0 );
  pp.t1  = sigproc.par_get('t1',  pars, 0 );
  pp.t2  = sigproc.par_get('t2',  pars, pp.t1 );
  pp.var       = sigproc.par_get('var',     pars, '' );
  pp.window    = sigproc.par_get('window',    pars, 0 );
  pp.remconst  = sigproc.par_get('remconst',  pars, 0 );
  do_plot      = sigproc.par_get('plot',      pars, 1 );
  list_file    = sigproc.par_get('list_file', pars, '' );
  % values in pp will go to cache

  if length(list_file)
    pdir = [ list_file '.trace'];
%    file_txt=[pdir '/' dstr '_' xfile pp.var '.txt'];
    file_png=[pdir '/' dstr '_' xfile pp.var '.png'];
    file_par=[pdir '/trace_data/' dstr '_' xfile pp.var '.mat'];
    if unix(['[ -s "' file_png '" -a -s "' file_par '" ]']) == 0
      load(file_par, '-mat', 'pars_cache', 'time', 'fre', 'amp');
      if isequal(pp, pars_cache)
        fprintf('skipping processed file: %s %s\n', dstr, xfile); 
        res=[time; fre; amp];
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/trace_data/']);
    % save original parameters before modifications
    pars_cache=pp;
  end

  % read signal:
  [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile);

  if pp.t2<=pp.t1; pp.t2=tx(end); end
  if pp.window<=0; pp.window=floor(length(tx)/100); end

  % cut signal t1 .. t2
  tx=tx-pp.t0;
  ii=find(tx>=pp.t1 & tx<=pp.t2);
  tx=tx(ii); xx=xx(ii);

%  % find starting frequency if it is not set yet
%  if pp.f0==0
%    [p2f,p2a]=sigproc.find_2peaks_st(tx,xx,pp.window, 10, pp.fmin, pp.fmax);
%    pp.f0=p2f(1);
%    [~, im] =max(p2a);
%    pp.f0=p2f(im);
%  end

  % find starting frequency if it is not set yet. 
  if pp.f0==0
    [f, a] = sigproc.fft(tx(1:pp.window), xx(1:pp.window),pp.fmin,pp.fmax);
    [ma, mi] = max(smooth(abs(a),10));
    pp.f0 = f(mi);
  end

  % print parameters
  fprintf('tracing parameters:\n');
  fprintf('  f0=%f df=%f t1=%f t2=%f window=%f remconst=%i plot=%i\n',...
     pp.f0, pp.df, pp.t1, pp.t2, pp.window, pp.remconst, do_plot);

  [T, F, A] = sigproc.fft_sl(tx, xx, pp.window, floor(pp.window/10),...
                             pp.fmin, pp.fmax);

  if pp.remconst==1;  A=sigproc.fft_sl_remconst(A, 20); end

  [fre,amp,dis] = sigproc.fft_sl_trace(F,A,pp.f0,pp.df,'max', pp.dfs);
  time=T; % this is needed for correct plot3d (time will be cutted)

  % cut data up to the end of traced signal
  i0=find(dis==0);
  if (length(i0)>0 && i0(1)>1)
    fre=fre(1:i0(1));
    amp=amp(1:i0(1));
    time=time(1:i0(1));
  end

  % plot data
  if do_plot
    ff=find_figure('aft'); clf; hold on;
    title([dstr ' ' xfile], 'Interpreter', 'none');
    subplot(2,2,1); hold on; title('amp - time');
      plot(time, amp, 'ob-', 'MarkerSize',2);
    subplot(2,2,2); hold on; title('freq - time');
      plot(time, fre, 'or-', 'MarkerSize',2);
    subplot(2,2,3); hold on; title('amp^2 - freq');
      plot(fre, amp.^2, 'om-', 'MarkerSize',2);
    subplot(2,2,4); hold on; title('3d');
      sigproc.plot_3di(T, F, abs(A), 'sqrt fft3d');
      plot (fre, time, 'g-');
      plot (fre-pp.df, time, 'g-');
      plot (fre+pp.df, time, 'g-');
      xlim([min(fre)-pp.df, max(fre)+pp.df]);
      ylim([min(time), max(time)]);
  else ff=0
  end

  % put txt and png into cache in the list mode
  if length(list_file)
    % save text file with time-fre-amp 
%    fo=fopen(file_txt, 'w');
%    fprintf(fo, '# %s %s\n', dstr, xfile); 
%    fprintf(fo, '# done by sig_trace.m program\n'); 
%    fprintf(fo, '# parameters: %s\n', pars); 
%    fprintf(fo, '#%10s %10s %10s\n', 'time', 'freq', 'amp'); 
%    for j=1:length(time)
%      fprintf(fo, '%f %f %f\n', time(j), fre(j), amp(j)); 
%    end
%    fclose(fo);

    % save png file 
    if ff>0
      set(ff, 'PaperUnits', 'points', 'PaperSize', [1280,1024],...
              'PaperPosition', [0,0,1280,1024]);
      figure(ff);
      print('-dpng', '-r72', file_png);
      unix(['chmod 664 ' file_png]);
    end

    save(file_par, 'pars_cache', 'time', 'amp', 'fre');
  end

  res=[time; fre; amp];

end
