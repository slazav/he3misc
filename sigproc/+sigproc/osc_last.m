function  [dstr, xfile] = last_sig()
%[dstr, xfile] = last_sig()
%
% get current date and name of last signal of this date

%%%% this function is obsoleted by sig_read.m

  [~, dstr] = unix('date "+%Y%m%d" | tr -d "\n"');
  [~, ddir] = unix('date "+/rota/data/%Y/%m/%d/osc" | tr -d "\n"');
  [i, xfile] = unix(['ls -1 -c ', ddir, ' | head -1 | tr -d "\n"']);
  if (i~=0) || strcmp(xfile,'')
    error('There are no signals in %s\n', ddir);
    return;
  end
  fprintf('%s %s\n', dstr, xfile); 
end