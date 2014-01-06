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
  pp.smin   = sigproc.par_get('smin',  pars, 2 );
  pp.smax   = sigproc.par_get('smax',  pars, pp.df );
  pp.exc_df = sigproc.par_get('exc_df', pars, 1 ); % filtering range for excitation signal

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
      res=load(file_par, '-mat');
      if usecache || (isequal(pp, res.pars) && isequal(pp_read, res.pars_read))
        fprintf('skipping processed file: %s %s\n', dstr, xfile);
        res.cache     = 1;
        return;
      end
    end
    unix(['mkdir -p -m 775 -- ' pdir '/freq/']);
  end


  % save original parameters before modifications
  pars_cache=pp;
  pars_read_cache=pp_read;
  time=[]; fre=[];

  % read signal
  if strcmp(xfile,'test') % model signal
      npts=2400000;
      osc_time=50; %s
      amp_sig0=0.123;
      fre_sig0=1543; %Hz
      mod_dep1=4.3; %Hz
      mod_dep2=17; %Hz
      mod_fre0=15; %Hz
    tx=linspace(0,osc_time, npts);
    dt_osc=(tx(end)-tx(1))/(npts-1);
    xx=amp_sig0 * sin(... % freq = d phase/ dt
         2*pi*fre_sig0.*tx -...
         mod_dep1/mod_fre0*cos(2*pi*mod_fre0*tx) -...
         mod_dep2/2/mod_fre0*cos(4*pi*mod_fre0*tx));

      amp_exc0=0.245;
      fre_exc0=15;
    txe=tx; dt_osce=dt_osc;
    xxe=amp_exc0 * sin(2*pi*fre_exc0*txe);
    xfile=sprintf('test_exc%.2fg%.2fHz_', 0.01, fre_exc0);
  else
    [dstr xfile ~] = sigproc.sig_ch_name(dstr, xfile);
    [tx,  xx, dt_osc]   = sigproc.sig_read(dstr, xfile, [pars ' chan=2']);
    [txe, xxe, dt_osce] = sigproc.sig_read(dstr, xfile, [pars ' chan=1']);
  end

  % print parameters
  fprintf('parameters: %s\n', pars);

  % set default window
  if pp.window<=0; pp.window=floor(length(tx)/100); end

  % do sliding fft
  if pp.step<=0; pp.step=round(pp.window/10); end
  [time,F,A] = sigproc.fft_sl(tx, xx, pp.window, pp.step, pp.fmin, pp.fmax);
  A=abs(A);

  % multiple signals
  if length(others); error('averaging is not supported\n'); end

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
  xxf = 2*ifft(X);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find exc freq and gain from the signal name
  exc_gain=0; exc_fre0=0;
  r = get_par(xfile, 'exc([0-9\.]*)g([0-9\.]*)Hz', 0);
  if r~=0; exc_gain=r(1); exc_fre0=r(2);
  else error('No exc<...>g<...>Hz text in the name!')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % fft of excitation signal
  [specfe, specae] = sigproc.fft(txe,xxe);
  ii=find(specfe>pp.smin & specfe<pp.smax);
  specfe=specfe(ii); specae=abs(specae(ii));

  % find maximum near exc_fre0 and get real
  % real frequency exc_fre, peak amplitude exc_amp and integral exc_int
  fmine=exc_fre0 - pp.exc_df;
  fmaxe=exc_fre0 + pp.exc_df;
  if fmine<pp.smin; fmine=pp.smin; end
  if fmaxe<pp.smin; fmaxe=pp.smin+pp.exc_df; end
  ii=find(specfe>fmine & specfe<fmaxe);
  [exc_amp,i]=max(specae(ii));
  exc_fre=specfe(ii(i));  % measured excitation freq
  exc_int=sqrt(sum(specae(ii).^2));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % filter exc signal
  fmine=exc_fre - pp.exc_df;
  fmaxe=exc_fre + pp.exc_df;
  if fmine<pp.smin; fmine=pp.smin; end
  if fmaxe<pp.smin; fmaxe=pp.smin+pp.exc_df; end
  ii_exc=find(specfe>fmine & specfe<fmaxe);

  N = length(txe);
  ddt = (txe(N)-txe(1))/(N-1);
  ddf = 1/ddt/N; % frequency step in fft
  X = fft(xxe);

  % apply filter
  ii = [round(fmine/ddf):round(fmaxe/ddf)];
  X1 = X(ii);
  X = zeros(1,N);
  X(ii) = X1;

  % reconstruct signal
  xxfe = 2*ifft(X);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find frequency from signal zeros
  tc = rel2f.find_zeros(tx,real(xxf));
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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % integrals of 1F, 2F, 4F, 6F peaks
  ii_int1=find(abs(specf-exc_fre)<pp.exc_df);
  sig_int1=sqrt(sum(speca(ii_int1).^2));

  ii_int2=find(abs(specf-2*exc_fre)<pp.exc_df);
  sig_int2=sqrt(sum(speca(ii_int2).^2));

  ii_int4=find(abs(specf-4*exc_fre)<pp.exc_df);
  sig_int4=sqrt(sum(speca(ii_int4).^2));

  ii_int6=find(abs(specf-6*exc_fre)<pp.exc_df);
  sig_int6=sqrt(sum(speca(ii_int6).^2));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % lock_in: multiply freq1 by xxfe
  sig1 = real(xxfe/exc_int);
  sig2 = 2*(real(xxfe/exc_int).^2 - mean(real(xxfe/exc_int).^2));
  % 2F signal is shifted by 1/4 period!

  mf=mean(freq1);
  freq2a = interp1(time1, freq1,txe) - mf;
  freq2b = interp1(time1 + 1/4/exc_fre, freq1,txe) - mf;
  freq2c = interp1(time1 + 1/8/exc_fre, freq1,txe) - mf;

  iia=find(~isnan(freq2a));
  iib=find(~isnan(freq2b));
  iic=find(~isnan(freq2c));
  sig_lia1 =  2*sum(freq2a(iia).*sig1(iia))/length(iia) +...
             2i*sum(freq2b(iib).*sig1(iib))/length(iib);
  sig_lia2 = 2i*sum(freq2a(iia).*sig2(iia))/length(iia) +...
              2*sum(freq2c(iic).*sig2(iic))/length(iic);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Simple lock-in: multiply by constructed signal with exc_fre freq
  % UPD. Bad idea: peak can be wider then 1/T, and we want integral but not maxinum
%  sig_slia1 = 2*abs(sum((freq1-mf).* exp(1i*2*pi*exc_fre*time1)))/length(time1);
%  sig_slia2 = 2*abs(sum((freq1-mf).* exp(1i*2*pi*2*exc_fre*time1)))/length(time1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot data
  if do_plot

    if length(list_file); t = 'sig_freq';
    else t = ['sig_freq: ' dstr ' - ' xfile pp.var];  end
    ff=find_figure(t); clf;
    tt=[min(time) max(time)];
    subplot(3,2,1); hold on;
      title([dstr ' ' xfile ' ' pp.var], 'Interpreter', 'none');
      ii = find(F > fmin & F < fmax);
      sigproc.plot_3di(F(ii), time, abs(A(ii,:))', 'sqrt');
      plot (tt, [fmin fmin], 'b-', 'MarkerSize',10);
      plot (tt, [fmax fmax], 'b-', 'MarkerSize',10);
      xlim(tt);
      ylim([fmin, fmax]);
    subplot(3,2,2); hold on;
      plot(time, freq, 'm-');
      xlim(tt);
      ylim([fmin, fmax]);
      ylabel('freq, Hz');
      xlabel('time, s');
    subplot(3,2,3:4); hold on;
      k=round(0.1*max(speca)/mean(speca));
      if k>4; plot(specf, k*speca, 'g-'); end
      plot(specf, speca, 'b-', 'linewidth', 2);
      plot(specf(ii_int1), speca(ii_int1), 'r-', 'linewidth', 2);
      plot(specf(ii_int2), speca(ii_int2), 'r-', 'linewidth', 2);
      %ke=max(speca)/max(specae)
      %plot(specfe, ke*specae, 'g-', 'linewidth', 2);
      ylim([0, max(speca)]);
      xlim([0, max(specf)]);
      ylabel('modulation amp, Hz');
      if k>4; legend(sprintf('x%d',k), 'x1'); end
      text(exc_fre+1, max(speca), {
          sprintf('Int: %.3e', sig_int1)
          sprintf('LI:  %.3e', abs(sig_lia1))
        },...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'top');
      text(2*exc_fre+1, max(speca), {
          sprintf('int: %.4e', sig_int2)
          sprintf('LI:  %.4e', abs(sig_lia2))
        },...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top');
      text(4*exc_fre+1, max(speca), {
          sprintf('int: %.4e', sig_int4)
%          sprintf('LI:  %.4e', abs(sig_lia4))
        },...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top');
      text(6*exc_fre+1, max(speca), {
          sprintf('int: %.4e', sig_int6)
%          sprintf('LI:  %.4e', abs(sig_lia4))
        },...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top');
    subplot(3,2,5:6); hold on;
      k=round(0.1*max(specae)/mean(specae));
      if k>4; plot(specfe, k*specae, 'g-'); end
      plot(specfe, specae,  'b-', 'linewidth', 2);
      plot(specfe(ii_exc), specae(ii_exc), 'r-', 'linewidth', 2);
      ylim([0, max(specae)]);
      xlim([0, max(specfe)]);
      xlabel('freq, Hz');
      ylabel('excitation amp');
      if k>4; legend(sprintf('x%d',k), 'x1'); end

      text(exc_fre+1, max(specae), {
          sprintf('Fset: %.2f', exc_fre0)
          sprintf('F:    %.2f', exc_fre)
          sprintf('A:    %.4e', exc_amp)
          sprintf('Int:  %.4e', exc_int)
        },...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top');

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


  res.exc_fre0  = exc_fre0; % set excitation frequency (from filename)
  res.exc_gain  = exc_gain; % set excitation gain (from filename)
  res.exc_fre   = exc_fre;  % measured excitation freq (peak maximum)
  res.exc_amp   = exc_amp;  % measured excitation amplitude (peak maximum)
  res.exc_int   = exc_int;  % integral of the excitation peak
  res.sig_int1  = sig_int1; % integral of signal frequncy 1F peak
  res.sig_int2  = sig_int2; % integral of signal frequncy 2F peak
  res.sig_lia1  = sig_lia1; % lock-in amplitude of 1F peak (complex)
  res.sig_lia2  = sig_lia2; % lock-in amplitude of 2F peak (complex)
  res.sig_slia1 = 0;        % simple lock-in amplitude of 2F peak
  res.sig_slia2 = 0;        % (leave this for compatibility with old data)
%  res.specf     = specf;
%  res.speca     = speca;
  res.pars      = pars_cache;
  res.pars_read = pars_read_cache;
  res.cache     = 0;

  % save results into png, save cache in the list mode
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
    % save mat file
    save(file_par, '-struct', 'res');
  end

end


  function res = get_par(fname, re, def)
    s= regexp(fname, re, 'tokens','once');
    if length(s)>0;
      for i=1:length(s)
        res(i)=str2num(s{i});
      end
    else res=def;
    end
  end
