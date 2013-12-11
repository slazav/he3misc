function data = sig_trace05(data, list_file)
  for i=1:length(data)

    % read + fft3di + trace parameters
    pars = sigproc2013.sig_trace05_pars(data{i}.pars);

    % other parameters are not saved in the cache
    retrace        = sigproc2013.par_get('retrace',     data{i}.pars, 0 ); % force processing
    do_plot        = sigproc2013.par_get('do_plot',     data{i}.pars, 1 ); % plot picture
    do_png         = sigproc2013.par_get('do_png',      data{i}.pars, 1 ); % create png files
    save_trace_txt = sigproc2013.par_get('save_trace_txt', data{i}.pars, '' ); % save result in a txt file
    save_trace_mat = sigproc2013.par_get('save_trace_mat', data{i}.pars, '' ); % save result in a mat file

    % retrace if there is no trace field
    % or parameters are different
    if ~isfield(data{i}, 'trace') ||...
       ~isfield(data{i}.trace, 'pars') || ~isequal(data{i}.trace.pars, pars)
       retrace=1;
    end

    % do retrace if we need png picture
    if length(list_file) && do_plot && do_png
      pdir = [ list_file '.cache/trace'];
      unix(['mkdir -p -m 775 -- ' pdir ]);
      file_png=[pdir '/' data{i}.alias '.png'];
      if unix(sprintf('[ -s "%s" ]', file_png)) ~= 0; retrace=1; end
    end

    if ~retrace; continue; end;

    % save original parameters before modifications
    data{i}.trace.pars = pars;

    % do fft
    [time, F, A, window, step] = sigproc2013.sig_fft3d(data{i}.dstr, data{i}.file, pars);

    % find starting frequency if it is not set yet
    if pars.f0<=0
      [~,mi] = max(abs(A(:,1)));
      pars.f0=F(mi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trace freq
    if pars.interactive
      find_figure('sig_trace - interactive'); clf; hold on;
      sigproc2013.plot_3di(F, time, abs(A)', 'sqrt');
      ylabel('frequency, Hz');
      xlabel('time, s');
      ylim([pars.fmin, pars.fmax]);
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
        if pars.ftracer==1
          if (j==1)  fre_ind=find(abs(F-pars.f0)<pars.df);
          else       fre_ind=find(abs(F-fre(j-1))<pars.df);
          end
          [~, mi] = max(abs(A(fre_ind,j)));
          mi=fre_ind(mi);
          fre(j)= F(mi);
        else
          if j==1;
            fre(j)=pars.f0; % current freq of a peak
            k0=find(F>=pars.f0,1); % k0 is an index of current peak position
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
              [fre(j),~,~] = sigproc2013.max_3pt(F((k0-1):(k0+1)), A((k0-1):(k0+1),j));
            end
          end
        end
        % check stopping condition
        fre_ind=find(abs(F-fre(j))<pars.df);
        N=length(fre_ind);
        win2=abs(linspace(1,-1,N)');
        win1=ones(N,1)-win2;
        ppp  = sum(abs(A(fre_ind,j)).^2 .* win1) / sum(win1) /...
               sum(abs(A(fre_ind,j)).^2 .* win2) * sum(win2);
        if (ppp<1 + pars.trace_th/100.0) break; end
      end
    end

    fend=mean(fre(max(1,j-10):max(1,j-1)));
    fre(j:length(time))=nan;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % integrate amp
    plot_fmin=max(min([pars.f0 fre])-pars.df, min(F));
    plot_fmax=min(max([pars.f0 fre])+pars.df, max(F));
    plot_ind=find(F > plot_fmin & F < plot_fmax);

    for j=1:length(time)
      %% set integrating range
      if isnan(fre(j)); f00=fend;
      else f00=fre(j);
      end

      if (pars.fixdf)
        sig_fmin(j)=plot_fmin;
        sig_fmax(j)=plot_fmax;
        sig_ind=plot_ind;
        ref_ind=1:length(F);
      else
        df=pars.df;
        if pars.autodf && j>1;
          % window length in seconds
          tw=(time(j)-time(j-1))*window/step;
          % linewidth from finite window, Hz
          df1 = 1/tw;
          % linewidth from frequency change
          ii=find(abs(time - time(j))<=tw/2);
          dff = abs(fre(ii(1))-fre(ii(end)));
          df2 = sqrt( 4/pi* dff/tw );
          df=8*max(df1,df2);
        end
        sig_fmin(j)=f00-df;
        sig_fmax(j)=f00+df;
        sig_ind=find(F>sig_fmin(j) & F<=sig_fmax(j));
        ref_ind=plot_ind;
      end
      %% integrate in sig_ind range and in the whole range
      amp2(j) = sum(abs(A(sig_ind,j)).^2);
      ampR(j) = sum(abs(A(ref_ind,j)).^2);
      noise2(j) = (ampR(j)-amp2(j))*length(sig_ind)/(length(ref_ind)-length(sig_ind));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot data
    if do_plot
      if length(list_file); t = 'sig_trace';
      else t = ['sig_trace: ' data{i}.alias ];  end
      ff=find_figure(t); clf;


      subplot(2,2,1); hold on;
        sigproc2013.plot_3di(F, time, abs(A)', 'sqrt');
        plot (time, sig_fmin, 'b-', 'MarkerSize',10);
        plot (time, sig_fmax, 'b-', 'MarkerSize',10);
        xlim([min(time), max(time)]);
        ylim([min(F), max(F)]);
        ylabel('frequency, Hz');
        xlabel('time, s');
      subplot(2,2,3); hold on;
        sigproc2013.plot_3di(F(plot_ind), time, abs(A(plot_ind,:))', 'sqrt');
        plot (time, fre, 'b-', 'MarkerSize', 10);
        plot (time, sig_fmin, 'b-', 'MarkerSize',10);
        plot (time, sig_fmax, 'b-', 'MarkerSize',10);
        xlim([min(time), max(time)]);
        ylim([plot_fmin, plot_fmax]);
        ylabel('frequency, Hz');
        xlabel('time, s');
      subplot(2,2,[2 4]); hold on;
        plot(time, log(amp2)/log(10), 'r-', 'MarkerSize', 10);
        plot(time, log(noise2)/log(10), 'b-', 'MarkerSize', 10);
        legend('log_{10}(amp^2)', 'log_{10}(noise^2)')
        xlim([min(time), max(time)]);
        xlabel('time, s');

      if iscell(data{i}.file); f=data{i}.file{1}; else f=data{i}.file; end
      f={data{i}.alias; f};
      text(0,1, f, 'units','normalized',...
         'Interpreter','none', 'HorizontalAlignment', 'center',...
         'VerticalAlignment', 'bottom');
    else ff=0
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save txt file
    if length(save_trace_txt)
      % save text file with time-fre-amp
      fo=fopen(save_trace_txt, 'w');
      fprintf(fo, '# done by sig_trace.m program\n');
      fprintf(fo, '# %s\n', data{i}.id);
      fprintf(fo, '# parameters: %s\n', pars);
      fprintf(fo, '#%10s %10s %10s %10s\n', 'time', 'freq', 'amp^2', 'noise^2');
      fprintf(fo, '%f %f %f %f\n', [time; fre; amp2; noise2]);
      fclose(fo);
    end
    % save mat file
    if length(save_trace_mat)
      id=data{i}.id;
      save(save_trace_mat, 'id', 'pars', 'time', 'fre', 'amp2', 'noise2');
    end

    % create png if needed 
    if length(list_file) && ff>0 && do_png
      set(ff, 'PaperUnits', 'points', 'PaperSize', [800,600],...
              'PaperPosition', [0,0,800,600]);
      print('-dpng', '-r72', file_png);
      unix(['chmod 664 ' file_png]);
    end

    data{i}.trace.time   = time;
    data{i}.trace.freq   = fre;
    data{i}.trace.amp2   = amp2;
    data{i}.trace.noise2 = noise2;
    data{i}.trace.window = window;
    data{i}.trace.step   = step;
  end
end
