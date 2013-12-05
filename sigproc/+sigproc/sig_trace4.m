function res = sig_trace4(dstr, xfile, pars)
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

  pp_read = sigproc.sig_read_pars(pars); % cache read parameters
  pp.fmin  = sigproc.par_get('fmin',  pars, 10 );
  pp.fmax  = sigproc.par_get('fmax',  pars, 3000 );
  pp.trfmin  = sigproc.par_get('trfmin',  pars, -1 );
  pp.trfmax  = sigproc.par_get('trfmax',  pars, -1 );
  pp.f0    = sigproc.par_get('f0',    pars, 0 );
  pp.df    = sigproc.par_get('df',    pars, 100 );
  pp.fixdf = sigproc.par_get('fixdf', pars, 0 );
  pp.ftracer  = sigproc.par_get('ftracer', pars, 2 );
  pp.fixnoise = sigproc.par_get('fixnoise', pars, 0 );
  pp.interactive = sigproc.par_get('interactive', pars, 0 );
  pp.fastfreq    = sigproc.par_get('fastfreq', pars, 0 );

  pp.fix_lock_in = sigproc.par_get('fix_lock_in', pars, 1 ); % fix amplitude according to lock-in curve
  pp.var         = sigproc.par_get('var',         pars, ''); % variant name (if there are several signals in one)
  pp.window      = sigproc.par_get('window',      pars, 0 );
  pp.step        = sigproc.par_get('step',        pars, 0 );

  % other parameterss are not saved in the cache
  do_plot        = sigproc.par_get('plot',        pars, 1 );
  list_file      = sigproc.par_get('list_file',   pars, '');
  refit          = sigproc.par_get('refit',       pars, 0 ); % don't use cache
  usecache       = sigproc.par_get('usecache',    pars, 0 ); % always use cache
  save_trace_txt = sigproc.par_get('save_trace_txt', pars, '' ); % save result in a txt file

  others={};
  if (iscell(xfile))
    others=xfile(2:end);
    xfile=xfile{1};
  end

  % get data from cache if possible
  if length(list_file)
    pdir = [ list_file '.cache'];
    file_png=[pdir '/trace/' dstr '_' xfile pp.var '.png'];
    file_par=[pdir '/trace/' dstr '_' xfile pp.var '.mat'];
    if ~refit && unix(['[ -f "' file_png '" -a -s "' file_par '" ]']) == 0
      load(file_par, '-mat', 'pars_cache', 'pars_read_cache',...
                             'time', 'fre', 'amp', 'wid', 'int2', 'noise2');
      if usecache || (isequal(pp, pars_cache) && isequal(pp_read, pars_read_cache))
        fprintf('skipping processed file: %s %s\n', dstr, xfile);
        res.time      = time;
        res.fre       = fre;
        res.amp       = amp;
        res.int2      = int2;
        res.noise2    = noise2;
        res.wid       = wid;
        res.pars      = pars_cache;
        res.pars_read = pars_read_cache;
        res.cache     = 1;
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/trace/']);
  end


  % save original parameters before modifications
  time=[]; amp=[]; fre=[]; wid=[]; int2 = []; noise2 = [];
  pars_cache=pp;
  pars_read_cache=pp_read;

  % read signal
  [tx, xx, dt_osc] = sigproc.sig_read(dstr, xfile, pars);

  % print parameters
  fprintf('tracing parameters: %s\n', pars);

  % set default window
  if pp.window<=0; pp.window=floor(length(tx)/100); end

  % do sliding fft
  if pp.step<=0; pp.step=round(pp.window/10); end
  [time,F,A] = sigproc.fft_sl(tx, xx, pp.window, pp.step, pp.fmin, pp.fmax);
  A=abs(A);

  % average multiple signals
  if length(others)
    A=A.^2;
    N=length(others);
    tl0=length(time);

    parfor i=1:N
      [tx1, xx1, dt_osc] = sigproc.osc_read(dstr, others{i}, pars);
      ii=find(tx1>=tx(1) & tx1<=tx(end));
      [time1,F1,A1{i}] = sigproc.fft_sl(tx1(ii), xx1(ii), pp.window, round(pp.window/10), pp.fmin, pp.fmax);
      tl(i) = length(time1);
    end

    tlm=min([tl tl0]);
    for i=1:N
      A = A(:,1:tlm) + abs(A1{i}(:,1:tlm)).^2;
    end
    time = time(1:min(tl));
    A=sqrt(A/(N+1));
  end

  % find starting frequency if it is not set yet
  if pp.f0<=0
    [~,mi] = max(abs(A(:,1)));
    pp.f0=F(mi);
  end

  % apply lock-in correction
  for j=1:length(time)
    if pp.fix_lock_in; A(:,j) = A(:,j)./sigproc.lock_in_gain(F)'; end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % trace freq
  if pp.interactive

    find_figure('sig_trace - interactive'); clf; hold on;
    sigproc.plot_3di(F, time, abs(A)', 'sqrt');
    ylabel('frequency, Hz');
    xlabel('time, s');
    ylim([pp.fmin, pp.fmax]);
    xlim([min(time), max(time)]);

    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')

    fii = []; tii=[];
    while 1
      [ti,fi,but] = ginput(1);
      if but~=1 break; end
      fii(end+1)=fi;
      tii(end+1)=ti;
      plot(tii,fii,'b*-');
    end
    [tii,nii] = sort(tii); ftt=fii(nii);
    fre = interp1(tii,fii, time);

  else
    for j=1:length(time)
      if pp.ftracer==1
        if (j==1)  fre_ind=find(abs(F-pp.f0)<pp.df);
        else       fre_ind=find(abs(F-fre(j-1))<pp.df);
        end
        [~, mi] = max(abs(A(fre_ind,j)));
        mi=fre_ind(mi);
        fre(j)= F(mi);
      else
        if j==1;
          fre(j)=pp.f0; % current freq of a peak
          k0=find(F>=pp.f0,1); % k0 is an index of current peak position
        else
          kmin=k0;
          kmax=k0;
          nmax=k0; % amp max on (kmin:kmax,j)
          while (1)
            a1=A(nmax,j);
            if kmin>1;         a2=A(kmin-1,j-1); else a2=a1-1; end
            if kmax<length(F); a3=A(kmax+1,j-1); else a3=a1-1; end
            [~,im] = max([a1 a2 a3]);
            % break if max in j larger then both sides in j-1
            if im==1; k0=nmax; fre(j)=F(k0); break; end
            % add largest side in j-1
            if im==2; kmin=kmin-1; add=kmin;
            else      kmax=kmax+1; add=kmax; end % add point
            % update max in j
            if A(nmax,j) < A(add,j); nmax=add; end
          end

          %additional fix
          if (k0>1 & k0<length(F))
            [fre(j),~,~] = sigproc.max_3pt(F((k0-1):(k0+1)), A((k0-1):(k0+1),j));
          end
%         fre_ind=find(abs(F-fre(j))<pp.df);
%         [~, k0] = max(abs(A(fre_ind,j)));
%         k0=fre_ind(k0);
%         fre(j)= F(k0);
        end
      end



      % check stopping condition
      fre_ind=find(abs(F-fre(j))<pp.df);
      N=length(fre_ind);

      win2=abs(linspace(1,-1,N)');
      win1=ones(N,1)-win2;

      ppp  = sum(abs(A(fre_ind,j)).^2 .* win1) / sum(win1) /...
             sum(abs(A(fre_ind,j)).^2 .* win2) * sum(win2);
      if (ppp<1.05) break; end
    end
  end

  fend=mean(fre(max(1,j-10):max(1,j-1)));
  fre(j:length(time))=nan;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % integrate amp
  if pp.trfmin < 0    
      plot_fmin=max(min([pp.f0 fre])-pp.df,pp.fmin);
  else
      plot_fmin = pp.trfmin;
  end
  
  if pp.trfmax < 0
      plot_fmax=min(max([pp.f0 fre])+pp.df,pp.fmax);
  else
      plot_fmax=pp.trfmax;
  end
  plot_ind=find(F > plot_fmin & F < plot_fmax);
  for j=1:length(time)

    %% set integrating range
    if isnan(fre(j)); f00=fend;
    else f00=fre(j);
    end

    if (pp.fixdf)
      sig_fmin(j)=plot_fmin;
      sig_fmax(j)=plot_fmax;
      sig_ind=plot_ind;
      ref_ind=1:length(F);
    else
      df=pp.df;
      if pp.fastfreq && j>1;
        df=max(df, abs(fre(j)-fre(j-1)) *pp.window/step / 2 );
      end
      sig_fmin(j)=f00-df;
      sig_fmax(j)=f00+df;
      sig_ind=find(F>sig_fmin(j) & F<=sig_fmax(j));
      ref_ind=plot_ind;
    end

    %% integrate in sig_ind range and in the whole range
    amp1(j) = sum(abs(A(sig_ind,j)).^2);
    amp2(j) = sum(abs(A(ref_ind,j)).^2);

    noise2(j) = (amp2(j)-amp1(j))*length(sig_ind)/(length(ref_ind)-length(sig_ind));
    if pp.fixnoise; amp(j) = amp1(j)-noise2(j);
    else amp(j)=amp1(j); end
  end
  if pp.fixnoise amp = amp + mean(noise2) - min(amp); end

  int2=amp1;

  amp = sqrt(amp);
  noise=sqrt(noise2);
  amp1=sqrt(amp1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot data
  if do_plot

    if length(list_file); t = 'sig_trace';
    else t = ['sig_trace: ' dstr ' - ' xfile];  end
    ff=find_figure(t); clf;

    subplot(2,2,1); hold on;
      title([dstr ' ' xfile ' ' pp.var], 'Interpreter', 'none');
      sigproc.plot_3di(F, time, abs(A)', 'sqrt');
%      plot (time, fre, 'b.-', 'MarkerSize',10);
      plot (time, sig_fmin, 'b-', 'MarkerSize',10);
      plot (time, sig_fmax, 'b-', 'MarkerSize',10);
      xlim([min(time), max(time)]);
      ylabel('frequency, Hz');
      xlabel('time, s');
      ylim([pp.fmin, pp.fmax]);
    subplot(2,2,3); hold on;
      sigproc.plot_3di(F(plot_ind), time, abs(A(plot_ind,:))', 'sqrt');
      plot (time, fre, 'b-', 'MarkerSize', 10);
      plot (time, sig_fmin, 'b-', 'MarkerSize',10);
      plot (time, sig_fmax, 'b-', 'MarkerSize',10);
      xlim([min(time), max(time)]);
      ylim([plot_fmin, plot_fmax]);
    subplot(2,2,[2 4]); hold on;
      plot(time, int2, 'm-', 'MarkerSize', 10);
      plot(time, noise2, 'b-', 'MarkerSize', 10);
      plot(time, int2-noise2, 'r.-', 'MarkerSize', 10);
      legend('amp^2', 'noise^2', 'amp^2-noise^2')
      xlim([min(time), max(time)]);
      ylabel('amplitude, V');
      xlabel('time, s');
  else ff=0
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % save ascii file
  if length(save_trace_txt)
    % save text file with time-fre-amp
    fo=fopen(save_trace_txt, 'w');
    fprintf(fo, '# %s %s\n', dstr, xfile);
    fprintf(fo, '# done by sig_trace.m program\n');
    fprintf(fo, '# %s %s\n', dstr, xfile);
    fprintf(fo, '# parameters: %s\n', pars);
    fprintf(fo, '#%10s %10s %10s\n', 'time', 'freq', 'amp');
    for j=1:length(time)
      fprintf(fo, '%f %f %f\n', time(j), fre(j), amp(j));
    end
    fclose(fo);
  end

  % put txt and png into cache in the list mode
  if length(list_file)

    % save png file.
    if ff>0
      set(ff, 'PaperUnits', 'points', 'PaperSize', [1280,1024],...
              'PaperPosition', [0,0,1280,1024]);
%      set(0,'CurrentFigure',ff);
      print('-dpng', '-r72', file_png);
      unix(['chmod 664 ' file_png]);
    else
      unix(['touch ' file_png]);
      unix(['chmod 664 ' file_png]);
    end

    save(file_par, 'pars_cache', 'pars_read_cache',...
                   'time', 'amp', 'fre', 'wid', 'int2', 'noise2');
  end

  res.time      = time;
  res.fre       = fre;
  res.amp       = amp;
  res.wid       = wid;
  res.int2      = int2;
  res.noise2    = noise2;
  res.pars      = pars_cache;
  res.pars_read = pars_read_cache;
  res.cache     = 0;
end

