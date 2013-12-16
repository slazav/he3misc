function sig_fit3d(varargin)
  %% Plot fft3d picture of the signal.
  %% In list operation data is not modified,
  %%   pictures are saved in cache folder.
  %% Can average signals.

  func=@sig_fft3d1;
  sigproc2013.sig_process_list(func, '', varargin{:});
end


function data=sig_fft3d1(data, list_file)

  % check png file 
  if length(list_file)
    pdir = [ list_file '.cache/fft3d'];
    unix(['mkdir -p -m 775 -- ' pdir ]);
    file_png=[pdir '/' data.alias '.png'];
    if unix(sprintf('[ -s "%s" ]', file_png)) == 0;
      fprintf('skipping processed file: %s\n', data.alias);
      return;
    end
  end

  [time, freq, amp] = sigproc2013.sig_fft3d(data.dstr, data.file, data.pars);

  if length(list_file); t = 'sig_fft3d';
  else t = ['sig_fft3d: ' data.alias];
  end

  ff=find_figure(t); clf; hold on;
  sigproc2013.plot_3di(time, freq, amp, data.pars);
  xlabel('freq, Hz');
  ylabel('time, s');
  if iscell(data.file); f=data.file{1}; else f=data.file; end
  f={data.alias; f};
  text(0.5,1, f, 'units','normalized',...
       'Interpreter','none', 'HorizontalAlignment', 'center',...
       'VerticalAlignment', 'bottom');

  % save png file 
  if length(list_file)
    figure(ff);
    set(ff, 'PaperUnits', 'points', 'PaperSize', [800 600],...
            'PaperPosition', [0 0 800 600]);
    print('-dpng', '-r72', file_png);
    unix(sprintf('chmod 664 "%s"', file_png));
  end
end
