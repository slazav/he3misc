function r=sig_fit(varargin)
%% Fit relaxation signals.
%
% sig_trace <file>  -- process signal list, write results to file.fit and file.fit_img
% sig_trace <date> <file>  -- proces one signal
% sig_trace                -- process last signal
% sig_trace <date> <file> <par> ...  -- use parameters (t1fit,t2fit,exp)

%  addpath('/rota/Analysis/PS/osc2011/'); % path with sigproc library

  func=@sigproc.sig_fit2;
  r=sigproc.foreach_sig(func, varargin{:});

%  % how to pass varargin directly to foreach_sig?? 
%  switch nargin
%    case 0; r=sigproc.foreach_sig(func);
%    case 1; r=sigproc.foreach_sig(func, varargin{1});
%    case 2; r=sigproc.foreach_sig(func, varargin{1}, varargin{2});
%    otherwise
%      pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
%        varargin(3:end), 'UniformOutput', false));
%      r=sigproc.foreach_sig(func, varargin{1}, varargin{2}, pars);
%  end


  if nargin==1 || (nargin>1 && length(varargin{2})==0)
    fn=varargin{1};
    fo=fopen([fn '.fit1'], 'w');
    n=length(r(1).file); % todo!!!

    fprintf( fo, '# This file was automatically generated from %s by sig_fit.m\n', fn); 
    fprintf(fo, ['%-8s %-' num2str(n) 's %-9s %-9s %-9s %-9s %-9s %-9s  %-9s %-9s %-9s  %s\n'], ... 
      '# date', 'file', 'A', 'aerr', 'tau,s', 'terr', 'b', 'berr', 'f0', 'df', 'tf', 'parameters'); 
    for i=1:length(r)
      fprintf(fo, ['%-8s %-' num2str(n) 's %9f %9f %9f %9f %9f %9f  %9.3f %9.3f %9.3f  %s\n'], ... 
         r(i).dstr, r(i).file, ...
         r(i).ampfit(1), r(i).ampfit(2), r(i).ampfit(3), r(i).ampfit(4), r(i).ampfit(5), r(i).ampfit(6),...
         r(i).frefit(1), r(i).frefit(2), r(i).frefit(3), r(i).pars );
    end
    fclose(fo);
  end
end

