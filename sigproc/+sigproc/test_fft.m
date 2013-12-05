function test_fft()

  tx=linspace(0,2.2,23102);
  xx=1.3*sin(2.0*pi*1001.2*tx);

  find_figure('test_fft'); clf; hold on;
  plot([995 1005],[0,0],'k-');

fprintf('sigproc.fft...\n');
  [fre1, amp1] = sigproc.fft(tx, xx, 0,2000);
  plot(fre1,real(amp1),'og');
  plot(fre1,imag(amp1),'ob');
  plot(fre1,abs(amp1),'or');

fprintf('sigproc.fft_prec...\n');
  [fre2, amp2] = sigproc.fft_prec(tx, xx, 995,1005, 0.05);
  plot(fre2,real(amp2),'g-');
  plot(fre2,imag(amp2),'b-');
  plot(fre2,abs(amp2),'r-');

fprintf('sigproc.peak_moments...\n');
  [a0, f0, d0] = sigproc.peak_moments(fre1, amp1, 1002, 100);
  plot(f0,a0,'k*');
  plot( f0-d0*[1,1], a0/2*[-1,1],'k-');
  plot( f0+d0*[1,1], a0/2*[-1,1],'k-');
  xlim([995,1005]);

end