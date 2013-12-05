function ilarm=get_ilarm(hmin, press)

%see /rota/Analysis/PS/osc2011/results_wz_20111207-20111213.dat
%see /rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/find_HL.m
%see /rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/find_HL_FS.m
  if press == 0.5
    Imin = [150 150 200 250 250 350 350 450 450 1500 1500 1250 1000];
    H_L  = [2637.563 2637.531 2638.705 2639.820 2639.861 2642.169 2642.195 ...
            2644.435 2644.468 2668.385 2668.354 2662.525 2656.808];
  elseif press == 0.0
    Imin = [ 250 250 950 1100 1400 450 165 250 250 300 350 400];
    H_L =  [2639 2639.05 2654.98 2658.5 2665.2 2643.52 2637 2638.955 ...
            2638.92 2639.998 2641.2 2642.365];
  else
    fprintf ('get_ilarm.m error: unknown pressure!\n');
    return;
  end
  P=polyfit(Imin,H_L,1);

  ilarm = polyval(P,hmin);
end