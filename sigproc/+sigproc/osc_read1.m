function [time, x, dt] = osc_read1(dstr, filename)
% Obsoleted by sig_read!

   % read raw oscilloscope signal (raw or .gz)

   fprintf('reading file: %s %s\n', dstr, filename);

    folder=['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/osc/'];
    infile=[folder filename];

    assert( unix(['[ -s "' infile '" ]']) == 0, '\nEmpty or missing: %s\n', infile);

    if regexp(filename, '.flac$') ~= 0
      f = [tempname '.wav'];
      unix(['flac -d  -s -o ' f ' -- ' infile ], '-echo');
      [data,fs,bits,opts] = wavread(f);
      % FLAC is stereo, the data is the right channel (channel No. 2)
      x = data(:,1)';
      dt = 1./fs;
      time = (0:length(data)-1) * dt;
      time = time - 0.1; % rough estimate of the burst is 100 ms after the start of sound card capture
      delete(f);
      return;
    elseif regexp(filename, '.wav$') ~= 0
      [data,fs,bits,opts] = wavread(filename);
      % WAV is stereo, the data is the right channel (channel No. 2) 
      x = data(:,1)';
      dt = 1./fs;
      time = (0:length(data)-1) * dt;
      time = time - 0.1; % rough estimate of the burst is 100 ms after the start of sound card capture
      return;
    elseif regexp(filename, '.gz$') ~= 0
      f = tempname;
      unix(['gunzip -c -- ' infile ' > ' f], '-echo');
      fin = fopen(f, 'r');
      delete(f);
    else
      fin = fopen(infile, 'r');
    end

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
%    if list(1) ~= 0
    jx = x<0;
    x(jx) = 256+x(jx);
      x=x-128;
%      x(find(x<-127)) = x(find(x<-127))+256;
%    end

    fclose(fin);
    x = x';
    x = x*dx+x0;
    %npts*dt;

    time = t0 + (0:dt:(length(x)-1)*dt);
%    time = 0:dt:(length(x)-1)*dt;

    fprintf('dt_osc: %g, points: %d\n', dt, length(time)); 
end