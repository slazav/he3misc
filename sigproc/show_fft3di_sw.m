function show_fft3di_sw(mode, dstr, xfile)

    addpath '/rota/Analysis/PS/osc2011/';

    colors='sqrt'; % color scale: normal | sqrt | log
    remconst=0;      % remove excitation signal: 0|1
    mode=str2num(mode);
%                    % 0 - no larmor calculations
                     % 1 - show larmor on the plot
                     % 2 - modify freq
                     % 3 - modify freq, calc sum
    width=2500;
    freq_corr = -400-(313.8-214.8);
    cut_stripe = 0;
%    cut_stripe = 0;
    cut_zero = 10; % Hz
    out_dir = '.';
    Fli = 827000; % Lock-in int freq

    if nargin == 1; [dstr, xfile] = sigproc.osc_last(); end

    fig_title = [dstr ' ' xfile];
    %fig_title = ['test'];
    [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile);

    % do fft
    window=round(0.8/dt_osc);
    step=round(0.1/dt_osc);
    [time, freq, amp] = sigproc.fft_sl(tx, xx, window, step, 10,width);
    % remove excitation signal
    if remconst==1; amp=sigproc.fft_sl_remconst(amp, 200); end
    if cut_zero>0 ; amp(find(freq<cut_zero),:) = 0; end

    if mode > 0
      fprintf('get larmor freq...\n');
      % get signal parameters
      [Time Iset Uexc Sens Omega Fork1 Hmin] = sigproc.osc_pars(xfile);
      % get NMR table
      [~, NMRt, ~, Iset, Imeas, ~, ~, ~, ~, ~, ~] =...
        sigproc.cwnmr_get(dstr, Time+time(1), Time+time(end));
      NMRt=NMRt-Time;
%Iset=Imeas;
      Iset = interp1q(NMRt',Iset',time'); % interpolate to time array
      % get f0
      [~, ~, f0] = sigproc.get_q(dstr,'A');
      FC1 = 49.04;  % freq correction for wrong osc/lock-in settings

      Ilarm = sigproc.get_ilarm(Hmin, 0.5); %get I_L from I_L vs Hmin data
  %    Ilarm=max(Iset); %sets the max current of the sweep to be I_L %BAD

      f_larm = f0*(1-(Ilarm-Iset).*Iset/(Ilarm^2)); % "real" larmor freq
      f_lc = FC1 + Fli - f_larm + freq_corr;  % larmor freq in our signal coords

      if (cut_stripe>0)
        jj=find(abs(freq-Fli+f0+FC1) > cut_stripe);
        amp(jj,:) = 0;
      end

    end
%%
    if mode > 1 % modify frequency
      fprintf('modify freq...\n');
      for i=1:length(time)
        ib=find(freq >=f_lc(i)-freq(end),1);
        ie=find(freq >=f_lc(i),1);
        if length(ie)==0; ie=length(freq); end
        if length(ib)==0; ib=1; end
        isrc=ib:ie;
        idst=ie:-1:ib;
        an = zeros(length(amp(:,i)),1);
        an(idst)=amp(isrc, i);
        amp(:, i) = an';
      end
    end
%%
    if mode > 2 % plot sum
      fprintf('plot sum...\n');
      find_figure(['sum:' fig_title]); clf; hold on;
      samp=sum(amp');
      samp=samp-smooth(freq, samp, 0.1, 'rloess')';
      plot(freq, samp, 'r-');

      %% save to file % tmp
      fo=fopen([out_dir '/' dstr '_' xfile '.sp'], 'w');
      fprintf(fo, '# %s %s\n', dstr, xfile);
      fprintf(fo, '# Time=%f Hmin=%f Omega=%f\n', Time, Hmin, Omega);
      for i=1:length(freq); fprintf(fo, '%f %f\n', freq(i), samp(i)); end
      fclose(fo);

      % find peaks
      [p_f, p_a] = rel2f.find_peaks(freq,samp);

      % select highest peaks
      [p_a, ii] = sort(p_a, 'descend');
      p_a=p_a(1:5);
      p_f=p_f(ii(1:5));
      % sort by freq
      [p_f, ii] = sort(p_f);
      p_a=p_a(ii);
      plot(p_f, p_a, 'm*');
      fprintf('f(0,0): %f\n', p_f(1));
    end

  if mode < 4
    find_figure(['3d: ' fig_title]); clf; hold on; title(fig_title);
    sigproc.plot_3di(time, freq, amp, 'scale=0 cbar=100');
    if mode==1; plot(f_lc, time, 'c-','LineWidth',2); end
  end


%    unix(['mkdir -p -- ', dstr]);
%    fbase=[dstr, '/' xfile];

%    hgsave([fbase '.3di.fig']);
%    print('-dpng', '-r150', [fbase, '.3di.png']);

%    lab=[fbase ...
%      '\nSIG:  dt_osc: ' num2str(dt_osc*1e6) ' mks, points: ' num2str(length(tx)) ', time: ' num2str(dt_osc*length(tx)) ' s' ...
%      '  FFT:  window: ' num2str(window*dt_osc) ' s, step: ' num2str(step*dt_osc) ' s'...
%      '\nSWEEP: ' num2str(sweep) 'mA/s' ...
%    ];
%    unix(['convert ', ...
%         ' label:"', lab, '" -gravity Center ',...
%         fbase, '.3di.png -append PNG8:', fbase, '.3di.png']);

end

