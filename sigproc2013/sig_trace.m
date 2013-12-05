function sig_trace(varargin)
%% Trace relaxation signals.
  func=@sigproc2013.sig_trace05;
  readonly=0;
  sigproc2013.sig_process_list(func, readonly, varargin{:});
end
