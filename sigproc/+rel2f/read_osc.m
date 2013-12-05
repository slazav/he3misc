function [time, x, dt] = read_osc(dstr, filename)

    folder=['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/osc/'];
    infile=[folder filename];

    if regexp(filename, '.gz$') ~= 0
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
    dx = list(8);
    x0 = list(9);
    for i=2:3
            fgetl(fin);
    end
    x = fscanf(fin, '%f');
    x = x';
    x = x*dx+x0;
    %npts*dt;
    time = t0 + (0:dt:(npts-1)*dt);
    fclose(fin);

