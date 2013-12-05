function sig_fit(varargin)
  func=@sig_fft_one;
  r=sigproc.foreach_sig(func, varargin{:});
end


function res=sig_fft_one(dstr, xfile, pars)

  [tx, xx, ~] = sigproc.sig_read(dstr, xfile, pars);
  list_file  = sigproc.par_get('list_file',  pars, '');
  var        = sigproc.par_get('var',        pars, '');
  fmin       = sigproc.par_get('minf',       pars, 10);
  fmin       = sigproc.par_get('fmin',       pars, fmin);
  fmax       = sigproc.par_get('maxf',       pars, 3000);
  fmax       = sigproc.par_get('fmax',       pars, fmax);

  [fre,amp] = sigproc.fft(tx, xx,fmin,fmax);

  if length(list_file); t = 'sig_fft';
  else t = ['sig_fft: ' dstr ' - ' xfile];
  end

  ff=find_figure(t); clf; hold on;
  ht=title([dstr ' - ' xfile]);
  set(ht,'Interpreter','none');

  plot(fre,abs(amp),'-r'); hold on;

  xlabel('freq, Hz');
  xlim([fre(1) fre(end)]);

  % save png file 
  if length(list_file)
    pdir = [ list_file '.cache/fft'];
    unix(['mkdir -p -m 775 -- ' pdir ]);
    file_png=[pdir '/' dstr '_' xfile var '.png'];
    set(ff, 'PaperUnits', 'points', 'PaperSize', [800 600],...
            'PaperPosition', [0 0 800 600]);
    figure(ff);
    print('-dpng', '-r72', file_png);
    unix(['chmod 664 ' file_png]);
  end
  res=ff;
end
