function sig_raw(varargin)
%% Plot raw signals.

  func=@sig_raw_one;

  % how to pass varargin directly to foreach_sig??
  switch nargin
    case 0; sigproc.foreach_sig(func);
    case 1; sigproc.foreach_sig(func, varargin{1});
    case 2; sigproc.foreach_sig(func, varargin{1}, varargin{2});
    otherwise
      pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
        varargin(3:end), 'UniformOutput', false));
      sigproc.foreach_sig(func, varargin{1}, varargin{2}, pars);
  end

end

function res=sig_raw_one(dstr, xfile, pars)

  [tx, x, ~] = sigproc.sig_read(dstr, xfile, pars);
  list_file  = sigproc.par_get('list_file',  pars, '');
  var        = sigproc.par_get('var',        pars, '');

  if length(list_file); t = 'sig_raw';
  else t = ['sig_raw: ' dstr ' - ' xfile];
  end
  ff=find_figure(t); clf; hold on;
  ht=title([dstr ' - ' xfile]);
  set(ht,'Interpreter','none');

  plot(tx,x);
  xlabel('time, s');
  xlim([tx(1) tx(end)]);

  % save png file
  if length(list_file)
    pdir = [ list_file '.cache/raw'];
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
