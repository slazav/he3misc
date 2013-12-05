function sig_process_list(func, readonly, dstr, xfile, varargin)
%% Run func for one signal or for signal list.
%% Use cached data from a .mat file

% sig_process_list <func> <listfile>           -- process signal list (dstr file pars)
% sig_process_list <func> <listfile> '' <pars> -- set additional parameters
% sig_process_list <func> <date> <file>        -- process one signal
% sig_process_list <func> <date> <file> <pars> -- process one signal with parameters
% sig_process_list <func>                      -- process last signal
%
% func is a function for datafile processing:
%   data = func(data, list_dile)
% data is a cell array of structures with following fields:
%   file   -- filename string or cell of strings
%   dstr   -- date string ('20131203') or cell of strings
%   folder -- folder ('/rota/data/2013/12/03') or cell of strings
%   id     -- signal id string (dstr + filename + number of averages + variant)
%   alias  -- alias string (to be written in plot headers etc)
%   avrg   -- number ov averages
%   pars   -- parameter string
%
% Parameters:
%   avrg -- Use averaging (0|1, default: 0).
%           All signals with same parameters (except timestamps)
%           will be averaged.
%   var  -- Signal variant (default '').
%           Use this parameter to process one signal multiple times.
%   skip -- Skip signal (0|1, default: 0).
%
% Special entries in the list
%  "common pars <parameters>" -- set parameters
%  "stop reading"

% 1. build a list of files and parameters
% 2. convert filenames, build ids and aliases
% 3. (if !readonly) read old datafile and transfer entries to the new list
% 4. apply func to each entry of the list
% 5. (if !readonly) save datafile

  data={};
  list_file='';
  data_file='';
  cmdline_pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
    varargin, 'UniformOutput', false));

  %%% 1. build a list of files and parameters
  if nargin == 3 || (nargin >= 4 && length(dstr)~=0 && length(xfile)==0)
    list_file=dstr;
    % list files contain three columns: dstr, filename, parameters
    [ndir, fn, pars] = textread(list_file, '%s %s %[^\n]\n',...
                       'commentstyle','shell');
    common_pars='';
    i=1;
    while i<=length(fn)
      d=char(ndir(i));
      f=char(fn(i));
      p=char(pars(i));

      % Read "common pars ..." line.
      if strcmp([d ' ' f], 'common pars')==1
        common_pars=[' ' p];
        i=i+1;
        continue;
      end

      if strcmp([d ' ' f], 'stop reading')==1
        break;
      end

      % parameters for current signal
      signal_pars=[cmdline_pars ' ' common_pars ' ' p];

      % Averaging!
      avrgn = 1;
      avrg  = sigproc2013.par_get('avrg', signal_pars, 0 );
      if avrg
        j=i+1;
        while j<=length(fn) % read next signals
          d1=char(ndir(j));
          f1=char(fn(j));
          % check that signal names after the timestamps are the same:
          t1 = regexp(f,  '(\_.*)', 'tokens', 'once');
          t2 = regexp(f1, '(\_.*)', 'tokens', 'once');
          if (~length(t1) || ~length(t2) || ~strcmp(t1{1},t2{1}))
            break;
          end
          j=j+1;
        end
        f=fn(i:j-1); % now f is a cell with many filenames
        d=ndir(i:j-1);
        avrgn=j-i;
        if (length(j)==0); break; end % ??
      end

      % Skipping signal
      skip  = sigproc2013.par_get('skip', signal_pars, 0 );
      if skip
        i=i+1;
        continue;
      end

      % Fill datalist
      last=length(data)+1;
      data{last}.file = f;
      data{last}.dstr = d;
      data{last}.pars = signal_pars;
      data{last}.avrg = avrgn;

      if avrg; i=j; else i=i+1; end
    end
  else  % one signal processing:
    if nargin == 1;
      dstr='last';
      xfile='last';
    end
    last=length(data)+1;
    data{last}.file = xfile;
    data{last}.dstr = dstr;
    data{last}.pars = cmdline_pars;
    data{last}.avrg = 1;
  end


  %%% 2. convert filenames, build ids and aliases
  data = sig_convert_filenames(data);

  %%% 3. read old datafile and transfer entries to the new list
  if ~readonly && length(list_file)
    data_file = [list_file '.mat'];
    if unix(sprintf('[ -s %s ]', data_file))==0
      old = load(data_file, 'data');
      for i=1:length(data)
        for j=1:length(old.data)
          if strcmp(data{i}.id, old.data{j}.id)
            old.data{j}.pars  = data{i}.pars; % update parameters
            old.data{j}.file  = data{i}.file;
            old.data{j}.dstr  = data{i}.dstr;
            old.data{j}.alias = data{i}.alias;
            data{i}=old.data{j};
          end
        end
      end
    end
  end

  %%% 4. apply func to each entry of the list
  data = func(data, list_file);

  %%% 5. save datafile
  if ~readonly && length(data_file)
     timestamp=now;
     save(data_file, 'data', 'list_file', 'timestamp');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = sig_convert_filenames(data)
  for i=1:length(data)
    % multiple files
    if iscell(data{i}.file)
      if ~iscell(data{i}.dstr) ||...
         length(data{i}.file) ~= length(data{i}.dstr)
        error('File and dstr cells have different size.\n');
      end
      for j=1:length(data{i}.file) % convert names
        [data{i}.dstr{j}, data{i}.file{j} data{i}.folder{j}] =...
          sigproc2013.sig_ch_name(data{i}.dstr{j}, data{i}.file{j});
      end
      % id is created from the first file:
      data{i}.id = [data{i}.dstr{1} ' ' data{i}.file{1}];
    else % single file
      % convert name:
      [data{i}.dstr, data{i}.file data{i}.folder] =...
        sigproc2013.sig_ch_name(data{i}.dstr, data{i}.file);
      data{i}.id = [data{i}.dstr ' ' data{i}.file];
    end
    % add number of averages to the id
    if data{i}.avrg>1
      data{i}.id = [ data{i}.id ' avrg:' num2str(data{i}.avrg)];
    end
    % add variant to the id
    var = sigproc2013.par_get('var', data{i}.pars, '');
    if length(var)
      data{i}.id = [ data{i}.id ' var:' var];
    end
    data{i}.alias=regexprep(data{i}.id,...
      '^([0-9]+ [0-9]+)[^ ]*', '$1');
    data{i}.alias=regexprep(data{i}.alias,...
      '[ :]', '_');
  end
end
