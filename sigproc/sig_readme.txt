    Signal processing

======================================================================
    Programs
sig_raw    - plot raw signal
sig_fft    - plot fft
sig_fft3d  - plot sliding fft
sig_update - update signal list
sig_trace  -
sig_fit    -
sig_forks  -

======================================================================
    Reading signals

[time, x, dt] = sigproc.sig_read(dstr, filename, pars)


dstr  -- date in YYYYMMDD format
         or 'last', or 'last-N'
xfile -- filename
         or 'last' or 'last-N' (signals are sorted by name)
         or time in [0-9]+ format

Signals are read from /rota/data/YYYY/MM/DD/osc/*.{osc,osc.gz,flac}.

    Parameters:

name    default
---------------------------------------------
auto0   0       t0 autodetection
chan    2       channel (for soundcard signals)
t0      0       time shift
t1      -inf    time range - begin
t2      +inf    time range - end
model   0       use model signal


This function is used in many "high-level" programs,
and this parameters and rules for signal names are available
there. Example:
  sig_raw last-1 last t1=0 t2=10 auto0
  sig_raw last 47562

======================================================================
    Processing lists and individual signals

Most of programs can process individual signals and signal lists.
This depends on arguments:

sig_raw <dstr> <xfile> <pars> ...
sig_raw <list>
sig_raw <list> '' <pars> ...

In list mode some data and plots can be created
in <list>.cache/ folder or in <list>.* files

======================================================================
    Plotting raw signal

sig_raw -- raw signal (amp vs time) is plotted.
In the list mode PNG picture is saved into
<list>.cache/raw/<dstr>_<xfile><var>.png
No data caching.

    Parameters:
-------------------------------------------------------------
var     ''       signal variant (used only in list mode to
                 process same file with different parameters)
+ read parameters (chan, t0,t1,t2, auto0, model)

======================================================================
    Plotting fft of a signal

sig_fft -- fourie spectrum of a signal (|amp|^2 vs freq) is plotted.
In the list mode PNG picture is saved into
<list>.cache/fft/<dstr>_<xfile><var>.png
No data caching.

    Parameters:
-------------------------------------------------------------
fmin    10       frequency range, min
fmax    3000     frequency range, max
var     ''       signal variant
+ read parameters (chan, t0,t1,t2, auto0, model)

======================================================================
    Plotting sliding fft of a signal

sig_fft3d

    Parameters:
-------------------------------------------------------------

fmin    10         frequency range, min
fmax    3000       frequency range, max
window  round(length(tx)/100) -- fft window
step    round(window/10)      -- sliding fft step
log     0          plot log(amp)
sqrt    1          plot sqrt(amp)
scale   0          maximal value (0 - auto)
colors  'copper2'  blue/copper/copper2/copper3
var     ''         signal variant
+ read parameters (chan, t0,t1,t2, auto0, model)

======================================================================


