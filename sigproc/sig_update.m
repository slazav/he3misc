function sig_update(file, varargin)
  % add signals to the list


  % join all additional arguments to pars
  pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
    varargin, 'UniformOutput', false));

  vars  = sigproc.par_get('vars',  pars, 0 );
  dt    = sigproc.par_get('dt',    pars, 0 );

  % find last dstr+xfile in a file
  fo=fopen(file);
  dat1=''; xfile='';
  dat2='';
  while ~feof(fo)
    l=fgets(fo);
    s = regexp(l, '^\s*(20[0-9]{2})([0-9]{2})([0-9]{2})\s+([^\s]*)','tokens','once');
    if length(s) >3;
      dat1=[s{1} '/' s{2} '/' s{3}];
      dat2=[s{1} s{2} s{3}];
      xfile=s{4};
    end
  end
  fclose(fo);

  % get newer files from this date
  xdir=['/rota/data/' dat1 '/osc/'];

  % problem: how to determin dir if there is nothing in the list?
  if length(xfile)
    newer = [' -newer ' xdir '/' xfile];
  else
    newer = '';
  end

%% this commented line just adds all new files to the list
%  unix([ 'find ' xdir newer -type f -not -empty -printf "' dat2 ' %f\n" | sort >> ' file ], '-echo');

%% more complicated thing: take vars and dt parameters into account:
  [ret, list] = unix([ 'find ' xdir newer [' -type f -not -empty -printf ' ...
                      '"'] dat2 ' %f\n" | sort -n -k2' ]);
  %disp(list);
  %return
  l1=regexp(list, '\n', 'split');

  for i=1:length(l1)
    sigfile = l1{i};
    if length(sigfile)==0; continue; end
    fo = fopen(file, 'a+');
    if (vars<=0)
      fprintf(fo, '%s\n', sigfile);
    else
      fprintf(fo, '\n');
      t0=0;
      for j=1:vars
        fprintf(fo, '%s var=%d t0=%f\n', sigfile, j, t0);
        t0=t0+dt;
      end
    end
    fclose(fo);
  end
end
