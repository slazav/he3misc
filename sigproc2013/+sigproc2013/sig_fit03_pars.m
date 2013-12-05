function pars = sig_read_pars(parstr)
% Parse parameters for sig_read.
% Only parameters which affects result are here.

    pars = sigproc2013.sig_trace05_pars(parstr); % read+fft3d+trace parameters
    % our own parameters

    % global fitting ranges
    pars.t1fit   = sigproc2013.par_get('t1fit',   parstr, 0 );  % time range
    pars.t2fit   = sigproc2013.par_get('t2fit',   parstr, pars.t1fit);
    pars.minamp  = sigproc2013.par_get('minamp',  parstr, 0 ); % amplitude range
    pars.maxamp  = sigproc2013.par_get('maxamp',  parstr, +inf );

    % freq fitting
    pars.t1fre      = sigproc2013.par_get('t1fre', parstr, pars.t1fit ); % time range for freq fitting
    pars.t2fre      = sigproc2013.par_get('t2fre', parstr, pars.t2fit );
    pars.minamp_fre = sigproc2013.par_get('minamp_fre',  parstr, pars.minamp ); % amplitude range for freq fitting
    pars.maxamp_fre = sigproc2013.par_get('maxamp_fre',  parstr, pars.maxamp );

    % amp-freq fitting
    pars.t1af       = sigproc2013.par_get('t1af', parstr, pars.t1fit ); % time range for amp-freq fitting
    pars.t2af       = sigproc2013.par_get('t2af', parstr, pars.t1fit );
    pars.minamp_af  = sigproc2013.par_get('minamp_af', parstr, pars.minamp ); % amplitude range for amp-freq fitting
    pars.maxamp_af  = sigproc2013.par_get('maxamp_af', parstr, pars.maxamp );

    % shrink fitting range using frequency change
    % (use positive value if frequency goes down)
    pars.maxdf_amp  = sigproc2013.par_get('maxdf_amp',   parstr, 0 );
    pars.maxdf_fre  = sigproc2013.par_get('maxdf_fre',   parstr, 0 );
    pars.maxdf_af   = sigproc2013.par_get('maxdf_af',    parstr, 0 );

    % fitting functions
    pars.func_amp = sigproc2013.par_get('func_amp', parstr, 0 ); % 0 - without base
    pars.func_fre = sigproc2013.par_get('func_fre', parstr, 1 ); % 0 - mean fre, 1 - exp.fit
    pars.func_af  = sigproc2013.par_get('func_af',  parstr, 0 );

    pars.fixnoise = sigproc2013.par_get('fixnoise',  parstr, 0 );
end
