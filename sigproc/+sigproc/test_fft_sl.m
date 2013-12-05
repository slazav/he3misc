function test_fft()

  tx=linspace(0,4.1,8000000);
  xx=(2.2-tx).*sin(2.0*pi*(400+300*tx) .*tx) +...
     1.3*cos(2.0*pi*1001*tx);


fprintf('sigproc.fft_sliding...\n');
  [t, f, a] = sigproc.fft_sl(tx, xx, 80000, 80000);
  a1 = sigproc.fft_sl_remconst(a);

  find_figure('fft_sl 3D'); clf; hold on; title('fft_sl 3D');
  sigproc.plot_3di(t, f, abs(a),  'sqrt');

  find_figure('fft_sl 3D1'); clf; hold on; title('fft_sl 3D1');
  sigproc.plot_3di(t, f, abs(a1), 'sqrt');

  find_figure('test_fft_sl'); clf; hold on;
  for i=1:length(t)
    plot(f, abs(a(:,i))+t(i),'r-');
    plot(f, abs(a1(:,i))+t(i),'b-');
  end
end
