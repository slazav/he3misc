function sig_raw(varargin)
  %% Plot raw picture of the signal.
  %% In list operation data is not modified,
  %%   pictures are saved in cache folder.

  func=@sig_raw1;
  sigproc2013.sig_process_list(func, '', varargin{:});
end


function data=sig_raw1(data, list_file)

  % multiple files are not supported here
  if iscell(data.file)
    error('sig_fft: averaging is not supported');
  end

  % check png file
  if length(list_file)
    pdir = [ list_file '.cache/raw'];
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

  % plot title
  if length(list_file); t = 'sig_raw';
  else t = ['sig_raw: ' data.alias];
  end

  % plot data
  ff=find_figure(t); clf; hold on;
  plot(tx,xx,'-b'); hold on;
  if iscell(data.file); f=data.file{1}; else f=data.file; end
  f={data.alias; f};
  text(0.5,1, f, 'units','normalized',...
       'Interpreter','none', 'HorizontalAlignment', 'center',...
       'VerticalAlignment', 'bottom');
  xlabel('time, s');

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
