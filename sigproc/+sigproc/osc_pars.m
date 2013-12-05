function [Time, Iset, Uexc, Sens, Omega, Fork1, Hmin] = osc_pars(fname)
% Obsoleted by osc_pars1.m
% get parameters from signal filename

% <prefix>_<time>s_<Iset>mA_<Uexc>mV_li<Sens>mV_<Omega>rad_<Fork1>mHz_e<Hmin>mA_sig.osc.gz

% <time>s_<prefix>_<Iset>mA_<Uexc>mV_li<Sens>mV_<Omega>rad_<Fork1>mHz_e<Hmin>mA_sig.osc.gz

    % new names
      s = regexp(fname,...
        '^([0-9.]*)s_.*_([0-9.]*)mA_([0-9.]*)mV_li([0-9.]*)mV_([0-9.-]*)rad_([0-9.]*)mHz_e([0-9.]*)mA',...
        'tokens','once');

    % old names
    if (length(s)==0)
      s = regexp(fname,...
        '_([0-9.]*)s_([0-9.]*)mA_([0-9.]*)mV_li([0-9.]*)mV_([0-9.-]*)rad_([0-9.]*)mHz_e([0-9.]*)mA',...
        'tokens','once');
    end

    % very old names
    if (length(s)==0)
      s = regexp(fname,...
        '_([0-9.]*)s_([0-9.]*)mA_([0-9.]*)mV_[0-9.]*sdiv_li([0-9.]*)mV_([0-9.-]*)rad_([0-9.]*)mHz_e([0-9.]*)m?A',...
        'tokens','once');
    end

    % short names
    if (length(s)==0)
      s = regexp(fname,...
        '^([0-9.]*)s_',...
        'tokens','once');
      s{2}='0'; s{3}='0'; s{4}='0'; s{5}='0'; s{6}='0'; s{7}='0';
    end


    for i=1:7; if length(s) < i; s{i} = '0'; end; end

    Time = str2num(s{1});
    Iset = str2num(s{2});
    Uexc = str2num(s{3});
    Sens = str2num(s{4});
    Omega = str2num(s{5});
    Fork1 = str2num(s{6});
    Hmin  = str2num(s{7});
end
