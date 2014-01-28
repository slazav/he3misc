function res = get_forks(d, f, p)

    if (iscell(f))
      for i=1:length(f)
        t(i)=get_time(d, f{i}, p);
      end
      Time=mean(t);
      f=f{1};
    else
      Time=get_time(d, f, p);
    end

    %% parameter: nextday -- use next day data (for full sweep signals written after midnight) 
    s=regexp(p, '(nextday)', 'tokens','once');
    if length(s)>0
      fmt='yyyymmdd';
      tc = datenum(d, fmt);
      d  = datestr(addtodate(tc, 1, 'day'), fmt);
    end

    %% parameter: midnight_tdep -- use tdep file saved to previous date
    s=regexp(p, '(midnight_tdep)','tokens','once');
    if length(s)>0
        fmt='yyyymmdd';
        tc = datenum(d, fmt);
        d = datestr(addtodate(tc,-1,'day'), fmt);
        Time = Time + 24*3600; %time values in fork-file continue over midnigth
    end

    %% parameter: procfork -- use procforc2 program 
    s=regexp(p, '(procfork)', 'tokens','once');
    ss=regexp(p, '(procfork1)','tokens','once'); %use procfork program. time limits should be defined!!
    if length(s)>0
      f1p=0; f2p=0;
%      f1 = procfork2('rotafork1',d, Time);
        if length(ss)>0
            fwid_t1=sigproc.par_get('fwid_t1',p,0);
            fwid_t2=sigproc.par_get('fwid_t2',p,90000);
            f1 = procfork('rotafork1',d,fwid_t1,fwid_t2);
            f2 = procfork('rotafork2',d,fwid_t1,fwid_t2);
        else
            f1 = procfork2('rotafork1',d, Time);
            f2 = procfork2('rotafork2',d, Time);
        end
      if size(f1)~=0; f1p=interp1(f1.tdep.time, f1.tdep.width, Time)*1000; else f1p=0; end
      if size(f2)~=0; f2p=interp1(f2.tdep.time, f2.tdep.width, Time)*1000; else f2p=0; end
      f1l=f1p; f2l=f2p; err1=0; err2=0;
    else % normal mode. 
      [f1p f1l err1] = sigproc.get_fork(d, Time, 'rotafork1');
      [f2p f2l err2] = sigproc.get_fork(d, Time, 'rotafork2');
    end
    res.f1p=f1p;
    res.f1l=f1l;
    res.f2p=f2p;
    res.f2l=f2l;
    res.T=Time;
    res.err1=err1;
    res.err2=err2;
    res.file=f;
    res.dstr=d;
end

function Time=get_time(d, f, p)
    %% parameter: sig_time=52111  -- time in seconds 
    %% or sig_time=12:02:04 
    Time=0;
    s=regexp(p, 'sig_time=([\d.]+)', 'tokens','once');
    if length(s)>0; Time = str2num(s{1}); end
    s=regexp(p, 'sig_time=(\d+):(\d+):(\d+)', 'tokens','once');
    if length(s)==3; Time = 3600*str2num(s{1}) + 60*str2num(s{2}) + str2num(s{3}); end
    if Time==0; % or get from filename 
      [Time Iset Uexc Sens Omega Fork1 Hmin] = sigproc.osc_pars(f);
    end
end
