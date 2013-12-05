function show_raw( dstr, xfile, pars )
%shows the raw signal from oscilloscope
% Obsoleted by sig_raw.m

%chan = sigproc.par_get('chan',      pars, 2 );

%read_func=@sigproc.osc_read;
%if chan==1 read_func=@sigproc.osc_read1;


if nargin < 2; [dstr, xfile] = sigproc.osc_last(); end
[tx, x, ~] = sigproc.osc_read(dstr, xfile);
find_figure(['Raw signal: ' xfile]);
plot(tx,x);

end

