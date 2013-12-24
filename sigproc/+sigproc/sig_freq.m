function res = sig_freq(dstr, xfile, pars)
% Frocess signals with frequency modulaton.
% Parameters:
% - same as for sig_read, sig_trace

  pp_read = sigproc.sig_read_pars(pars); % cache read parameters
  pp.fmin  = sigproc.par_get('fmin',  pars, 10 );
  pp.fmax  = sigproc.par_get('fmax',  pars, 3000 );
  pp.fmin  = sigproc.par_get('minf',  pars, pp.fmin );
  pp.fmax  = sigproc.par_get('maxf',  pars, pp.fmax );
  pp.f0    = sigproc.par_get('f0',    pars, 0 );
  pp.df    = sigproc.par_get('df',    pars, 100 );
  pp.ftracer  = sigproc.par_get('ftracer', pars, 2 );

  pp.fix_lock_in = sigproc.par_get('fix_lock_in', pars, 1 ); % fix amplitude according to lock-in curve
  pp.var         = sigproc.par_get('var',         pars, ''); % variant name (if there are several signals in one)
  pp.window      = sigproc.par_get('window',      pars, 0 );
  pp.step        = sigproc.par_get('step',        pars, 0 );

  % spectrum range
  pp.smin  = sigproc.par_get('smin',  pars, 2 );
  pp.smax  = sigproc.par_get('smax',  pars, 120 );

  % integration ranges
  pp.ir12a  = sigproc.par_get('ir12a',  pars, 12 );
  pp.ir12b  = sigproc.par_get('ir12b',  pars, 13 );
  pp.ir25a  = sigproc.par_get('ir25a',  pars, 24 );
  pp.ir25b  = sigproc.par_get('ir25b',  pars, 26 );
  pp.ir50a  = sigproc.par_get('ir50a',  pars, 49 );
  pp.ir50b  = sigproc.par_get('ir50b',  pars, 51 );

  % other parameterss are not saved in the cache
  do_plot        = sigproc.par_get('plot',        pars, 1 );
  list_file      = sigproc.par_get('list_file',   pars, '');
  refit          = sigproc.par_get('refit',       pars, 0 ); % don't use cache
  usecache       = sigproc.par_get('usecache',    pars, 0 ); % always use cache
  save_freq_txt = sigproc.par_get('save_freq_txt', pars, '' ); % save result in a txt file

  others={};
  if (iscell(xfile))
    others=xfile(2:end);
    xfile=xfile{1};
  end

  % get data from cache if possible
  if length(list_file)
    pdir = [ list_file '.cache'];
    file_png=[pdir '/freq/' dstr '_' xfile pp.var '.png'];
    file_par=[pdir '/freq/' dstr '_' xfile pp.var '.mat'];
    if ~refit && unix(['[ -f "' file_png '" -a -s "' file_par '" ]']) == 0
      load(file_par, '-mat');
      if usecache || (isequal(pp, pars_cache) && isequal(pp_read, pars_read_cache))
        fprintf('skipping processed file: %s %s\n', dstr, xfile);
        res.specf     = specf;
        res.speca     = speca;
        res.a12       = a12;
        res.a25       = a25;
        res.a50       = a50;
        res.pars      = pars_cache;
        res.pars_read = pars_read_cache;
        res.cache     = 1;
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/freq/']);
  end


  % save original parameters before modifications
  time=[]; amp=[]; fre=[]; wid=[]; int2 = []; noise2 = [];
  pars_cache=pp;
  pars_read_cache=pp_read;

  % read signal
  [tx, xx, dt_osc] = sigproc.sig_read(dstr, xfile, pars);

  % print parameters
  fprintf('parameters: %s\n', pars);

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
  % do usual freq tracing
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

  fend=mean(fre(max(1,j-10):max(1,j-1)));
  fre(j:length(time))=nan;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % filter signal according to the found frequency range
  fmin=max(min([pp.f0 fre])-pp.df,pp.fmin);
  fmax=min(max([pp.f0 fre])+pp.df,pp.fmax);

  N = length(tx);
  ddt = (tx(N)-tx(1))/(N-1);
  ddf = 1/ddt/N; % frequency step in fft
  X = fft(xx);

  % apply filter
  ii = [round(fmin/ddf):round(fmax/ddf)];
  X1 = X(ii);
  X = zeros(1,N);
  %X(ii) = X1 .* blackman(length(X1))';
  X(ii) = X1;

  % reconstruct signal
  xxf = real(ifft(X));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find frequency from signal zeros
  tc = rel2f.find_zeros(tx,xxf);
  time1 = (tc(1:end-1)+tc(2:end))/2;
  freq1 = 0.5 ./ diff(tc);

  % longest range which is within fmin-fmax
  ii1 = find(~isnan(freq1) & freq1>=fmin & freq1<=fmax);
  vv = [0 cumsum(diff(ii1)~=1)];
  ii2 = ii1(vv==mode(vv));
  time1=time1(ii2); freq1=freq1(ii2);

  % interpolate on uniform time grid
  nfft = 2^(nextpow2(length(freq1))+1);
  time = linspace(time1(1),time1(end),nfft);
  freq = interp1(time1,freq1,time,'spline');

  % do FFT
  specf = 0.5/diff(time([1 2]))*linspace(0,1,nfft/2+1);
  speca = fft(freq,nfft)/nfft;
  speca = abs(2*speca(1:nfft/2+1)); %peak-amplitude spectrum, single-sided
  ii=find(specf>=pp.smin & specf<=pp.smax);
  specf=specf(ii); speca=speca(ii);

  % find peak integrals
  power=(speca/sqrt(2)).^2;
  ii12=find(specf>=pp.ir12a & specf <=pp.ir12b);
  ii25=find(specf>=pp.ir25a & specf <=pp.ir25b);
  ii50=find(specf>=pp.ir50a & specf <=pp.ir50b);
  a12 = sqrt(sum(power(ii12)));
  a25 = sqrt(sum(power(ii25)));
  a50 = sqrt(sum(power(ii50)));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot data
  if do_plot

    if length(list_file); t = 'sig_freq';
    else t = ['sig_freq: ' dstr ' - ' xfile];  end
    ff=find_figure(t); clf;
    tt=[min(time) max(time)];
    subplot(2,2,1); hold on;
      title([dstr ' ' xfile ' ' pp.var], 'Interpreter', 'none');
      ii = find(F > fmin & F < fmax);
      sigproc.plot_3di(F(ii), time, abs(A(ii,:))', 'sqrt');
      plot (tt, [fmin fmin], 'b-', 'MarkerSize',10);
      plot (tt, [fmax fmax], 'b-', 'MarkerSize',10);
      xlim(tt);
      ylim([fmin, fmax]);
    subplot(2,2,2); hold on;
      plot(time, freq, 'm-');
      xlim(tt);
      ylim([fmin, fmax]);
      ylabel('freq, Hz');
      xlabel('time, s');
    subplot(2,2,3:4); hold on;
      k=round(0.1*max(speca)/mean(speca));
      if k>4; plot(specf, k*speca, 'b-'); end
      plot(specf, speca, 'r-', 'linewidth', 2);
      ylim([0, max(speca)]);
      xlabel('freq, Hz');
      ylabel('modulation amp, Hz');
      if k>4; legend(sprintf('x%d',k), 'x1'); end
      text(10, max(speca), num2str(a12),...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom');
      text(25, max(speca), num2str(a25),...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom');
      text(50, max(speca), num2str(a50),...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom');

  else ff=0
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % save ascii file
  if length(save_freq_txt)
    % save text file with time-fre-amp
    fo=fopen(save_freq_txt, 'w');
    fprintf(fo, '# %s %s\n', dstr, xfile);
    fprintf(fo, '# done by sig_freq.m program\n');
    fprintf(fo, '# %s %s\n', dstr, xfile);
    fprintf(fo, '# parameters: %s\n', pars);
    fprintf(fo, '#%10s %10s %10s\n', 'specf', 'speca');
    for j=1:length(time)
      fprintf(fo, '%f %f\n', specf(j), speca(j));
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
                   'speca', 'specf', 'a12', 'a25', 'a50');
  end

  res.speca     = speca;
  res.specf     = specf;
  res.a12       = a12;
  res.a25       = a25;
  res.a50       = a50;
  res.pars      = pars_cache;
  res.pars_read = pars_read_cache;
  res.cache     = 0;
end

