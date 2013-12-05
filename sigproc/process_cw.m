function process_cw(dstr, t1, t2)
  addpath /rota/Analysis/PS/osc2011;

  chan='A';
  sdir=-1;
  t1=str2num(t1);
  t2=str2num(t2);
  % try to find additional information in CW.txt file
  [dstra, t1a, t2a, pressa, f0, hmina, rota, fwa, chana, sdira, ~] =...
    textread('cw_list.txt', '%s %f %f %f %f %f %f %f %s %f%[^\n]',...
      'commentstyle', 'shell');
  press=nan; hmin=nan; rot=nan; fw=nan;

  for i = 1:length(dstra)
    if strcmp(dstra{i},dstr) && t1a(i)==t1 && t2a(i)==t2
      press=pressa(i); hmin=hmina(i); rot=rota(i);
      fw=fwa(i); chan=chana{i}; sdir=sdira(i); break;
    end
  end

  pars=sprintf('chan=%s dir=%f', chan, sdir)

  % read and average spectra
  [I A N] = sigproc.cwnmr_avrg(dstr, t1,t2, pars);

  % plot spectra for processing
  sigproc.plot_spectrum(I,A, [ 'data/' dstr '-' num2str(t1) '-' num2str(t2)]);

end