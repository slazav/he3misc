function show_fft3di_list(varargin)

  addpath('/rota/Analysis/PS/osc2011/'); % path with sigproc library

  % how to pass varargin directly to foreach_sig?? 
  prog = @show_fft3di;
  switch nargin
    case 0; r=sigproc.foreach_sig(prog);
    case 1; r=sigproc.foreach_sig(prog, varargin{1});
    case 2; r=sigproc.foreach_sig(prog, varargin{1}, varargin{2});
    otherwise
      pars = cell2mat(cellfun(@(x) horzcat(x, ' '),...
        varargin(3:end), 'UniformOutput', false));
      r=sigproc.foreach_sig(prog, varargin{1}, varargin{2}, pars);
  end

end

