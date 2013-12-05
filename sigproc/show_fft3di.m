function h = show_fft3di(dstr, xfile, varargin)
%output is the handle for 3d-like image of the data

  if nargin < 1; dstr=''; end
  if nargin < 2; xfile=''; end

  pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
         varargin, 'UniformOutput', false));

  if length(dstr)==0 || length(xfile)==0; [dstr, xfile] = sigproc.osc_last(); end

  % get parameters:
  window = sigproc.par_get('window',  pars, 0 );
  step   = sigproc.par_get('step',  pars, 0 );
  minf   = sigproc.par_get('minf',  pars, 20 );
  maxf   = sigproc.par_get('maxf',  pars, 3000 );
  t0     = sigproc.par_get('t0',    pars, 0 );
  t1     = sigproc.par_get('t1',    pars, -inf );
  t2     = sigproc.par_get('t2',    pars, t1 );
  scale  = sigproc.par_get('scale', pars, 'sqrt' );
  remconst   = sigproc.par_get('remconst',  pars, 0 );
  list_file  = sigproc.par_get('list_file', pars, '' );
  chan       = sigproc.par_get('chan',      pars, 2 );

%  read_func=@sigproc.osc_read;
%  if chan==1 read_func=@sigproc.osc_read1; end

  fig_title = [dstr ' ' xfile ' ' chan];
%  [tx, xx, dt_osc] = read_func(dstr, xfile);

  [tx, xx, dt_osc] = sigproc.osc_read(dstr, xfile, pars);

  % set default t1, t2 and window 
  if t1<=tx(1); t1=tx(1); end
  if t2<=t1; t2=tx(end); end
  if window==0; window=floor(length(tx)/100); end
  if step==0; step=floor(window/10); end

  % subtract t0, cut to [t1 .. t2] 
  tx=tx-t0;
  ii=find(tx>=t1 & tx<=t2);
  tx=tx(ii); xx=xx(ii);


  [time, freq, amp] = sigproc.fft_sl(tx, xx, window, step, minf, maxf);
  xx = []; tx = [];
  amp=abs(amp);

  if remconst; amp=sigproc.fft_sl_remconst(amp, 20); end

  find_figure(['3d: ' fig_title]); clf; hold on; title(fig_title);
  h = sigproc.plot_3di(time, freq, amp, pars);

  if length(list_file)>0
    mydir = [list_file '.trace/3d'];
    unix(['mkdir -p -m 775 -- ' mydir]);
    fbase=[mydir '/' dstr '_' xfile];
    print('-dpng', '-r150', [fbase, '.png']);
%    hgsave([fbase, '.3d.fig']);
  end
end

