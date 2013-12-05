function res = osc_pars(fname)
% get parameters from signal filename

  function res = get_par(fname, re, def)
    s= regexp(fname, re, 'tokens','once');
    if length(s)>0;
      for i=1:length(s)
        res(i)=str2num(s{i});
      end
    else res=def;
    end
  end

  res.time = get_par(fname, '^([0-9\.]*)s_', 0);

  res.fl    = get_par(fname, '_li([0-9\.]*)Hz_', 0);
  res.fp    = get_par(fname, '_pi([0-9\.]*)Hz',  0);

  if res.fl==0; res.fl = get_par(fname, '_([0-9\.]*)Hz', 0); end
  if res.fp==0; res.fp = get_par(fname, '_([0-9\.]*)Hz', 0); end

  res.rotsync = get_par(fname, '_sy([0-9\.]*)s_',  0);
  res.ncyc    = get_par(fname, '([0-9\.]*)ncyc',   0);
  res.vpulse  = get_par(fname, 'Hz([0-9\.]*)V',    0);

  res.iset = get_par(fname, '_([0-9\.]*)mA_', 0);
  if res.iset==0; res.iset = get_par(fname, '_ia([0-9\.]*)mA_', 0); end

  res.omega = get_par(fname, '_([0-9\.]*)rad_', 0);
  if res.omega==0; res.omega = get_par(fname, '_om([0-9\.]*)rs_', 0); end

  res.hmin = get_par(fname, '_e([0-9\.]*)mA', 0);

  res.sens = get_par(fname, '_li([0-9\.]*)mV_', 0);
  res.fw = get_par(fname, '_([0-9\.]*)mHz', 0);

  res.dem = get_par(fname, '_demag([0-9\.]*)A', 0);
  
  
  r = get_par(fname, '_sound([0-9\.]*)V([0-9\.]*)Hz', 0);
  if r~=0
  res.snd_amp=r(1); res.snd_fre=r(2);
  end
end
