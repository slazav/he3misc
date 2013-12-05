function val = par_get(name, pars, def)
% Get parameter value.
% Look arguments for <name> and return <value> for "<name>=<value>" substring,
%  1 for "<name>", 0 for "no_<name>", <def> if nothing found.
% First found value is used. Multiple parameters in one argument must be
% separated by whitespaces. Output type (char/double) depends on <def> type!

  if length(pars)==0; val=def; return; end

  s=regexp(pars, ['(^|\s+)(no_)?' name '(=(\S+))?($|\s+)'], 'tokens','once');

  if length(s)~=4; val=def; return; end
  if strcmp(s{2}, 'no_'); val=0; return; end
  if strcmp(s{3}, ''); val=1; return; end
  val=s{3}(2:end);

  if strcmp(class(def), 'double'); val=str2num(val); end
end
