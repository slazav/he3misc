function sig_fit(varargin)
%% Fit relaxation signals.
%
% sig_trace <file>  -- process signal list, write results to file.fit and file.fit_img
% sig_trace <date> <file>  -- proces one signal
% sig_trace                -- process last signal
% sig_trace <date> <file> <par> ...  -- use parameters (t1fit,t2fit,exp)

  addpath('/rota/Analysis/PS/osc2011/'); % path with sigproc library

  func=@sigproc.get_forks;

  % how to pass varargin directly to foreach_sig?? 
  switch nargin
    case 0; r=sigproc.foreach_sig(func);
    case 1; r=sigproc.foreach_sig(func, varargin{1});
    case 2; r=sigproc.foreach_sig(func, varargin{1}, varargin{2});
    otherwise
      pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
        varargin(3:end), 'UniformOutput', false));
      r=sigproc.foreach_sig(func, varargin{1}, varargin{2}, pars);
  end

  if nargin==1 || (nargin>1 && length(varargin{2})==0)
    fn=varargin{1};
    fo=fopen([fn '.frk'], 'w');
    n=length(r(1).file); % todo!!!

    fprintf( fo, '# This file was automatically generated from %s by sig_forks.m\n', fn); 
    fprintf( fo, ['%-8s %-' num2str(n) 's %7s %7s  %7s %7s  %11s %7s %7s\n'], ... 
      '# date', 'file', 'f1p,mHz', 'f1l,mHz', 'f2p,mHz', 'f2l,mHz', 'time,s', 'err1', 'err2'); 

    for i=1:length(r)
      fprintf( fo, ['%-8s %-' num2str(n) 's %7.2f %7.2f  %7.2f %7.2f  %11.5f %7.4f %7.4f\n'], ... 
         r(i).dstr, r(i).file, r(i).f1p, r(i).f1l, r(i).f2p, r(i).f2l, r(i).T, r(i).err1, r(i).err2);

    end
    fclose(fo);
  end

end

