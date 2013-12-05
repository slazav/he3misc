function foreach_sig_test()
  dstr='20120913';
  file='s4a826700.0Hz4.0V5000cyc_67280.88654s_2633.20mA_50.0mV_li10mV_1.02rad_418.76mHz_e250.0mA_sig.osc.gz';

  % sigproc.foreach_sig(@f); % last

  fprintf('One signal with default parameters:\n');
  sigproc.foreach_sig(@f, dstr, file);

  fprintf('One signal with p1=p1 p2=2:\n');
  sigproc.foreach_sig(@f, dstr, file, 'p1=p1 p2=2');

  fprintf('The same:\n');
  sigproc.foreach_sig(@f, dstr, file, 'p1=p1', 'p2=2');

  fprintf('List:\n');
  sigproc.foreach_sig(@f, 'foreach_sig_test.txt');

  fprintf('List with overwrited p1=p1new:\n');
  sigproc.foreach_sig(@f, 'foreach_sig_test.txt', '', 'p1=new');

end


function res = f(dstr, file, pars)
  fprintf('  dstr: %s\n  file: %s\n', dstr, file);
  fprintf('  pars: %s\n', pars);
  p1 = sigproc.par_get('p1', pars, 'p1def'); % string parameter
  p2 = sigproc.par_get('p2', pars,  -1.0);   % numerical parameter
  fprintf('  p1: %s, p2: %f\n\n', p1, p2);
  res = [p1 length(p2)];
end