function sig_fft(varargin)
  %% Plot fft picture of the signal.
  %% In list operation data is not modified,
  %%   pictures are saved in cache folder.

  func=@sig_fft1;
  sigproc2013.sig_process_list(func, '', varargin{:});
end


function data=sig_fft1(data, list_file)

  % multiple files are not supported here
  if iscell(data.file)
    error('sig_fft: averaging is not supported');
  end

  % check png file
  if length(list_file)
    pdir = [ list_file '.cache/fft'];
    unix(['mkdir -p -m 775 -- ' pdir ]);
    file_png=[pdir '/' data.alias '.png'];
    if unix(sprintf('[ -s "%s" ]', file_png)) == 0;
      fprintf('skipping processed file: %s\n', data.alias); 
      return;
    end
  end

  % read file
  pars = data.pars;
  [tx, xx, ~] = sigproc2013.sig_read(...
    data.dstr, data.file, data.pars);

  fmin = sigproc2013.par_get('minf', pars, 10);
  fmin = sigproc2013.par_get('fmin', pars, fmin);
  fmax = sigproc2013.par_get('maxf', pars, 3000);
  fmax = sigproc2013.par_get('fmax', pars, fmax);

  [fre,amp] = sigproc2013.fft(tx,xx,fmin,fmax);

  % plot title
  if length(list_file); t = 'sig_fft';
  else t = ['sig_fft: ' data.alias];
  end

  % plot data
  ff=find_figure(t); clf; hold on;
  plot(fre,abs(amp),'-r'); hold on;
  if iscell(data.file); f=data.file{1}; else f=data.file; end
  f={data.alias; f};
  text(0.5,1, f, 'units','normalized',...
       'Interpreter','none', 'HorizontalAlignment', 'center',...
       'VerticalAlignment', 'bottom');
  xlabel('freq, Hz');
  xlim([fre(1) fre(end)]);

  % save png file 
  if length(list_file)
    figure(ff);
    set(ff, 'PaperUnits', 'points', 'PaperSize', [800 600],...
            'PaperPosition', [0 0 800 600]);
    print('-dpng', '-r72', file_png);
    unix(sprintf('chmod 664 "%s"', file_png));
  end
  res=ff;
end
