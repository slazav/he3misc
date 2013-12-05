function h = show_fft3di_date(dstr)

    sigproc.check_dstr(dstr);

    folder=['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/osc/'];
    out_dir=[ dstr, '_3d'];

    list=dir(folder);

    unix(['mkdir -p -m 775 -- ' out_dir ]);

    for i=1:length(list)
      f = list(i).name;

      if f(1)=='.'; continue; end

      png = [out_dir '/' f '.png'];

      if unix(['[ -s "' png '" ]']) == 0;
        fprintf('%s %s\n', dstr, f); 
        continue;
      end

      fig_title = [dstr ' ' f];
      [tx, xx, dt_osc] = sigproc.osc_read(dstr, f);

      window=floor(length(tx)/100)*2;
      step=floor(window/10);

      if (tx(end)-tx(1) > 2)
        minf = 100;
        maxf = 1600;
      else
        minf = 824000;
        maxf = 827500;
%        minf = 170000;
%        maxf = 180000;
      end

      [time, freq, amp] = sigproc.fft_sl(tx, xx, window, step, minf, maxf);
      amp=abs(amp);

%save('tmp.mat', 'time', 'freq','amp');

      find_figure('3d'); clf; hold on; title(fig_title);
      h = sigproc.plot_3di(time, freq, amp, 'sqrt');
      print('-dpng', '-r150', png);
%      hgsave(fig);
    end
end

