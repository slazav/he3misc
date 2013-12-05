function [time, fre, amp, window, step] = sig_fft3d(dstr, xfile, pars)
% Do sliding FFT of a number of signals and return averaged (mean square) amplitude.
% Dstr and xfile are strings or cell arrays with strings.
% Parallel calculations can be used
%
% -- slazav, dec 2013.

  if nargin<3; pars=''; end
  if ~isstruct(pars)
    pars = sigproc2013.sig_fft3d_pars(pars);
  end

  if ~iscell(xfile)
    [tx, xx, dt_osc] = sigproc2013.sig_read(dstr, xfile, pars);
    [time,fre,amp,window,step] = sigproc2013.fft_sl(tx, xx,...
      pars.window, pars.step, pars.fmin, pars.fmax);
    % apply lock-in correction
    amp=abs(amp);
    for j=1:length(time)
      if pars.fix_lock_in;
        amp(:,j) = amp(:,j)./sigproc2013.lock_in_gain(fre);
      end
    end
    return
  end

  N=length(xfile);
  t={}; f={}; a={}; tmin=[];
  parfor i=1:N
    [tx, xx, dt_osc] = sigproc2013.sig_read(dstr{i}, xfile{i}, pars);
    [t{i},f{i},a{i},window,step] = sigproc2013.fft_sl(tx, xx,...
      pars.window, pars.step, pars.fmin, pars.fmax);
    tmin(i)=t{i}(1);
  end

  % starting time, and number of points
  % for each signal
  t1=max(tmin); imin=[]; nmax=[];
  for i=1:N
    imin(i) = find(t{i}>=t1,1);
    if length(imin(i))
      nmax(i) = length(t{i}) - imin(i);
    else
      nmax(i) = 0;
    end
  end
  nt = min(nmax(find(nmax>0)));

%  if phase>0
%    for k=2:length(time);
%      dt=time(k)-time(k-1);
%      amp(:,k) = amp(:,k).*exp(-i*freq'/2.0/pi *dt);
%    end
%  end
%  if     phase==1 amp=real(amp);
%  elseif phase==2 amp=imag(amp);
%  end
%  amp=abs(amp);

  time = t{1}(imin(1):imin(1)+nt);
  fre  = f{1};
  amp  = zeros(length(fre), nt+1);
  for i=1:N
    amp = amp + abs(a{i}(:, imin(i):imin(i)+nt)).^2;
  end
  amp = sqrt(amp/N);
  % apply lock-in correction 
  for j=1:length(time)
    if pars.fix_lock_in;
      amp(:,j) = amp(:,j)./sigproc2013.lock_in_gain(fre);
    end
  end
end
