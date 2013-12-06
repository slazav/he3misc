function pars = sig_fft3d_pars(parstr)
% Parse parameters for sig_fft3d + sig_read.
% Only parameters which affects result are here.

  pars = sigproc2013.sig_read_pars(parstr); % read parameters
  % our own parameters
  pars.fmin   = sigproc2013.par_get('minf',       parstr, 10);
  pars.fmin   = sigproc2013.par_get('fmin',       parstr, pars.fmin);
  pars.fmax   = sigproc2013.par_get('maxf',       parstr, 3000);
  pars.fmax   = sigproc2013.par_get('fmax',       parstr, pars.fmax);
  pars.window = sigproc2013.par_get('window',     parstr, -1);
  pars.step   = sigproc2013.par_get('step',       parstr, -1);
  pars.phase  = sigproc2013.par_get('phase',      parstr, 0);  % 0 - show amp, 1 - show real part. 
  pars.fix_lock_in = sigproc2013.par_get('fix_lock_in', parstr, 1);

end
