function [time, x, dt] = sig_read(dstr, xfile, pars)
  % read raw oscilloscope/soundcard signal (raw or .gz or flac.gz)

  if nargin<3; pars=''; end
  pp = sigproc.sig_read_pars(pars);

  [dstr, xfile, folder] = sigproc.sig_ch_name(dstr, xfile);
  fprintf('reading file: %s %s\n', dstr, xfile);

  infile=[folder xfile];
  assert( unix(['[ -s "' infile '" ]']) == 0, '\nEmpty or missing: %s\n', infile);

  if pp.model==1 %% model signal
    t00 = 3.01;
    f00 = 1512.0;
    df00 = 100.0;
    tf00 = 4.5;
    a00 = 1.1;
    a00n = 0;
    time=linspace(0,15, 3000000);
    dt=(time(end)-time(1))/length(time);
      t0m = t00*(1+0.5*sin(2*pi*500*time));
    amp00 = a00*exp(-time./t0m);
    fre00 = f00 + df00*exp(-time/tf00);
    ph00 = 2*pi*( time*f00 - df00*tf00*exp(-time/tf00) );
    noise = random('normal', 0, a00n, size(time));
    x=amp00.*sin(ph00) + noise;
    fmt=0;
  elseif pp.model==2 %% model signal
    f00 = 1512.0;
    f01 = 100;
    a00 = 1.1;
    time=linspace(0,15, 3000000);
    dt=(time(end)-time(1))/length(time);
    ph00 = 2*pi*( time*f00 - f01*time.^2);
    x=a00.*sin(ph00);
    fmt=0;
  elseif regexp(xfile, '.flac$') ~= 0
    f = [tempname '.wav'];
    unix(['flac -d  -s -o ' f ' -- ' infile ], '-echo');
    [data,fs,bits,opts] = wavread(f);
    % FLAC is stereo, the data is the right channel (channel No. 2)
    x = data(:,pp.chan)';
    dt = 1./fs;
    time = (0:length(data)-1) * dt;
    delete(f);
    fmt=1;
  elseif regexp(xfile, '.wav$') ~= 0
    [data,fs,bits,opts] = wavread(xfile);
    % WAV is stereo, the data is the right channel (channel No. 2) 
    x = data(:,pp.chan)';
    dt = 1./fs;
    time = (0:length(data)-1) * dt;
    fmt=1;
  elseif regexp(xfile, '.gz$') ~= 0
    f = tempname;
    unix(['gunzip -c -- ' infile ' > ' f], '-echo');
    fin = fopen(f, 'r');
    delete(f);
    fmt=2;
  else
    fin = fopen(infile, 'r');
    fmt=2;
  end

  if fmt==2 % osc signal
    head = fgetl(fin);
    list = sscanf(head,'%d,%d,%d,%d,%e,%e,%d,%e,%e,%d');
    npts = list(3);
    dt = list(5);
    t0 = list(6);
    ts = list(7);
    dx = list(8);
    x0 = list(9);
    xs = list(10);
    for i=2:3; fgetl(fin); end
    x = fscanf(fin, '%d');

    % 28/05/2013 - check header before converting binary data between
    % signed and unsigned int.
    if list(1) == 0
      jx = x<0;
      x(jx) = 256+x(jx);
      x=x-128;
    end
    fclose(fin);
    x = x';
    x = x*dx+x0;
    time = t0 + (0:dt:(length(x)-1)*dt);
  end

  %% find pulse in channel pp.auto0 (for soundcard) or in chan 2 (osc)
  if pp.auto0
    if fmt==1
      ii=find(time>=pp.auto0st);
      am=max(abs(data(ii,pp.auto0)))
      im=find(abs(data(ii,pp.auto0)) > am*pp.auto0th, 1)
    else
      [~, im] = max(abs(x));
    end
    time = time-time(im);
  end

  %% apply t0,t1,t2 parameters
  if pp.t0 time=time-pp.t0; end
  ii=find(time>=pp.t1 & time<=pp.t2);
  x=x(ii); time=time(ii);

  fprintf('dt_osc: %g, points: %d\n', dt, length(time));
end
