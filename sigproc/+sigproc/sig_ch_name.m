function [dstr xfile folder] = sig_ch_name(dstr, xfile)
%% convert special dir and filenames

%% Empty string or word ``last'' in {\tt dstr} is converted to current date
%% string (YYYYMMDD). ``last-{\it N}'' is converted to the date which was
%% {\it N} days ago.

  if length(dstr)==0;  dstr='last'; end

  s = regexp(dstr, '^last(-[0-9]*)?$', 'tokens','once');
  if (length(s))
    if length(s{1}) d=[' -d "' s{1} ' days"'];
    else d=''; end
    [~, dstr] = unix(['date "+%Y%m%d" ' d ' | tr -d "\n"']);
  end

%% {\tt folder} is set to {\tt /rota/data/YYYY/MM/DD/osc/'
%% according to converted dstr.

  folder=['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/osc/'];

%% Empty string or word ``last'' in {\tt xfile} is converted to the last
%% file in the {\tt folder}. ``last-{\it N}'' is converted to
%% {\it (N+1)}-th file from the end. If {\tt xfile} if a number, it is
%% converted to the filename which starts from this number
%% (it should be a timestamp).

  if length(xfile)==0; xfile='last'; end
  s = regexp(xfile, '^last(-[0-9]*)?$', 'tokens','once');
  if (length(s))
    if length(s{1})
      N = str2num(s{1}(2:end))+1;
      flt = ['sed -n -e "' num2str(N) 'p"'];
    else flt='head -1'; end
    [i, xfile] = unix(['ls -1 -c ', folder, ' | ' flt ' | tr -d "\n"']);
  end

  s = regexp(xfile, '^([0-9]+)$', 'tokens','once');
  if (length(s))
    [i, xfile] = unix(['ls -1 -c ', folder, ' | grep -m 1 "^' s{1} '" | tr -d "\n"']);
  end
end
