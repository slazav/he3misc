function pars = sig_trace05_pars(parstr)

    pars = sigproc2013.sig_fft3d_pars(parstr); % read+fft3d parameters
    % our own parameters
    pars.f0          = sigproc2013.par_get('f0',       parstr, -1 ); % starting freq (default: max amp at a first row)
    pars.df          = sigproc2013.par_get('df',       parstr, 30 ); % frequency window
    pars.trace_th    = sigproc2013.par_get('trace_th', parstr, 20 );  % trace threshold, percent
    pars.interactive = sigproc2013.par_get('interactive', parstr, 0 );
    pars.ftracer     = sigproc2013.par_get('ftracer', parstr, 2 ); % 1-simple (maximum inside df range), 2-smart
    pars.fixdf       = sigproc2013.par_get('fixdf',   parstr, 0 );  % integrate amplitude in a constant range, calculate noise in the whole range!
    pars.autodf      = sigproc2013.par_get('autodf',  parstr, 0 );
end
