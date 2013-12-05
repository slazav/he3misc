function h = plot_3di_v2(time, freq, amp, scale, fig_title,f0)
% plot_3d(time, freq, amp, scale, title)
%
% plot 3d-like picture of a signal
% scale can be 'auto' or real number
%
% slazav, feb 2012

    %% plot 3d-like frequency picture 

    amp=abs(amp);
    amin=min(min(amp));

    if strcmp(scale, 'log')      amp=log(amp-amin);
    elseif strcmp(scale, 'sqrt') amp=sqrt(amp-amin);
    end

    dat=[
      0.0  1 1 1
      0.2  0 0 1
      0.4  0 1 1
      0.6  0 1 0
      0.8  1 1 0
      1.0  1 0 0 ];

    dat_copper=[
      0.0  0 0 0
      0.5  1 0 0
      1.0  1 1 0 ];

    dat_copper1=[
      0.0    0 0 0
      1/3.0  1 0 0
      2/3.0  1 1 0
      1.0    1 1 1 ];

    dat_copper2=[
      0.0  1 1 1
      0.4  0 0 1
      0.6  1 0 1
      0.8  1 0 0
      1.0  1 1 0 ];

    amin=min(min(amp));
    amax=max(max(amp));
    
       
    
    [freqp,timep]=meshgrid(f0-freq,time);
    amp(amp<10^-2*max(amp(:)))=0;
    
     s.time=timep;
   s.amp=amp;
   s.freq=freqp;
   
   save('ampl_vs_time_plot_sautti.mat','s');
   
    h(1) = surf(freqp, timep,...
       (amp'-amin),'LineStyle','none');
   
   colormap hot  
   shading flat
   %shading interp
   
   ht=title(fig_title);
   set(ht,'Interpreter','none');
   
   
end

function col = rb(a, dat)
  v=dat(:,1); r=dat(:,2); g=dat(:,3); b=dat(:,4);
  col=zeros(size(a,2),size(a,1),3);
  col(:,:,1) = interp1(v, r, a', 'linear', 0); % R
  col(:,:,2) = interp1(v, g, a', 'linear', 0); % G
  col(:,:,3) = interp1(v, b, a', 'linear', 0); % B
  col(find(col>1))=1;
  col(find(col<0))=0;
end

