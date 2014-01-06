function res=sig_trace(varargin)
%% Trace relaxation signals.
%  addpath('/rota/Analysis/PS/osc2011/'); % path with sigproc library

  func=@sigproc.sig_freq;
  res=sigproc.foreach_sig(func, varargin{:})

end
