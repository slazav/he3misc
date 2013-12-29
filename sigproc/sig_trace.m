function r=sig_trace(varargin)
%% Trace relaxation signals.
%
% sig_trace <file>           -- process signal list, write results to file.trace dir
% sig_trace <file> '' <pars> -- process signal list, write results to file.trace dir
% sig_trace <date> <file>    -- proces one signal
% sig_trace                  -- process last signal
% sig_trace <date> <file> <pars> -- use parameters (t1,t2,f0,df,remconst)

%  addpath('/rota/Analysis/PS/osc2011/'); % path with sigproc library

  func=@sigproc.sig_trace4;
  r=sigproc.foreach_sig(func, varargin{:})

end
