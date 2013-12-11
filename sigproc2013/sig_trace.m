function sig_trace(varargin)
%% Trace relaxation signals.
  func=@sigproc2013.sig_trace05;
  sigproc2013.sig_process_list(func, 'readonly=0', varargin{:});
end
