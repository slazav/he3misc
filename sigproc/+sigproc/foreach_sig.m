function res = foreach_sig(func, dstr, xfile, varargin)
%% Run func for each signal.
%
% foreach_sig <func> <file>           -- process signal list (dstr file pars)
% foreach_sig <func> <file> '' <pars> -- set additional parameters
% foreach_sig <func> <date> <file>    -- proces one signal
% foreach_sig <func>                  -- process last signal
% foreach_sig <func> <date> <file> <pars> -- use parameters
%
% func is a function for single signal processing:
%   res = func(dstr, xfile, pars)
% Res can be numbers, vectors, of matrixes with constant height.
% If list file is used then its name is passed as additional
% parameter list_file=<file>, and number of the signal in list_num=<n>.


  % signal list processing
  if nargin == 2 || (nargin >= 3 && length(dstr)~=0 && length(xfile)==0)
    pars_add = cell2mat(cellfun(@(x) horzcat(x, ' '),...
      varargin, 'UniformOutput', false));


    [ndir, fn, pars] = textread(dstr,'%s %s %[^\n]\n',...
                       'commentstyle','shell');
    res=[];
    common_pars='';
    i=1;
    while i<=length(fn)
      d=char(ndir(i));
      f=char(fn(i));
      p=char(pars(i));

      if strcmp([d ' ' f], 'common pars')==1
        common_pars=[' ' p];
        i=i+1;
        continue;
      end

      avrg  = sigproc.par_get('avrg',   [pars_add common_pars], 0 );

      if avrg
        j=i+1;
        while j<=length(fn)
          d1=char(ndir(j));
          f1=char(fn(j));
          t1 = regexp(f,  '(\_.*)', 'tokens', 'once');
          t2 = regexp(f1, '(\_.*)', 'tokens', 'once');
%          if (~strcmp(f(12:end),f1(12:end)))
          if (~length(t1) || ~length(t2) || ~strcmp(t1{1},t2{1}))
            break;
          end
          j=j+1;
        end
        f=fn(i:j-1);
        p = [ p ' var=_avrg' num2str(j-i) ];
        if (length(j)==0); break; end
      end

      u = func(d, f, [pars_add ' ' p ' list_file=' dstr ' list_num=' num2str(i) common_pars]);
      res = horzcat(res, u);
      if avrg; i=j; else i=i+1; end
    end
  else  % one signal processing:
    if nargin == 1; [dstr, xfile] = sigproc.osc_last(); end % last signal

    pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
      varargin, 'UniformOutput', false));
    res=func(dstr, xfile, pars);
  end

end

