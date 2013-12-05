% this is a modification of par_get function from sigproc library
% which works for both matlab and octave
function val = par_get(name, pars, def)
% Get parameter value.
% Look arguments for <name> and return <value> for "<name>=<value>" substring,.
%  1 for "<name>", 0 for "no_<name>", <def> if nothing found..
% First found value is used. Multiple parameters in one argument must be.
% separated by whitespaces. Output type (char/double) depends on <def> type!.

  if length(pars)==0; val=def; return; end

  s=regexp(pars, ['(^|\s+)(no_)?' name '(=(\S+))?($|\s+)'], 'tokens','once');
  if length(s)==0; val=def; return; end

  for i=1:length(s)
    if strcmp(s{i}, 'no_'); val=0; return; end
    if length(s{i})==0; continue; end

    if strcmp(s{i}(1), '=');
      val=s{i}(2:end);
      if strcmp(class(def), 'double'); val=str2num(val); end
      return;
    end
  end
  val=1;
  return;
end
