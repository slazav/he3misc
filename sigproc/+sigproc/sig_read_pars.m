function pars = sig_read_pars(parstr)
  % Parse parameters for sig_read.
  % Only parameters which affects result are here.

  pars.chan    = sigproc.par_get('chan',    parstr, 2);    % channel
  pars.auto0   = sigproc.par_get('auto0',   parstr, 0);    % pulse autodetection channel
  pars.auto0th = sigproc.par_get('auto0th', parstr, 0.9);  % pulse autodetection threshold
  pars.auto0st = sigproc.par_get('auto0st', parstr, -inf); % starting time of pulse autodetection
  pars.t0      = sigproc.par_get('t0',      parstr, 0);    % time shift
  pars.t1      = sigproc.par_get('t1',      parstr, -inf); % time range - begin
  pars.t2      = sigproc.par_get('t2',      parstr, +inf); % time range - end
  pars.model   = sigproc.par_get('model',   parstr, 0);    % use model signal
end
