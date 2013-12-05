function h = plot_3di(time, freq, amp, pars)
% plot_3d(time, freq, amp, pars)
%
% plot 3d-like picture of a signal
%
% slazav, feb 2012

    %% plot 3d-like frequency picture 

    if nargin < 4; pars=''; end

    logscale  = sigproc.par_get('log',    pars, 0 );
    sqrtscale = sigproc.par_get('sqrt',   pars, 1 );
    cutamp    = sigproc.par_get('cutamp', pars, -1 );
    colors    = sigproc.par_get('colors', pars, 'copper2' );
    cbar      = sigproc.par_get('cbar',   pars, 5 );
    swapxy    = sigproc.par_get('swapxy', pars, 0 );
    flarm     = sigproc.par_get('flarm',  pars, 0 );

    if swapxy>0; amp=amp'; end
    amp=abs(amp);
    if (cutamp>0); amp(find(amp>cutamp)) = cutamp; end


    amin0=min(min(amp));
    amax0=max(max(amp));

    if cbar
      data = linspace(amin0,amax0,length(amp(1,:)));
      amp(end-cbar+1:end,:)=repmat(data,cbar,1);
      amp(end-cbar,:)=amin0;
    end

    if flarm; freq=flarm-freq; end

    if logscale;  amp=log(amp-amin0); end
    if sqrtscale; amp=sqrt(amp-amin0); end

    if strcmp(colors, 'blue');
      dat=[
        0.0  1 1 1
        0.2  0 0 1
        0.4  0 1 1
        0.6  0 1 0
        0.8  1 1 0
        1.0  1 0 0 ];
    elseif strcmp(colors, 'copper');
      dat=[
        0.0  0 0 0
        0.5  1 0 0
        1.0  1 1 0 ];
    elseif strcmp(colors, 'copper2');
      dat=[
        0.0    0 0 0
        1/3.0  1 0 0
        2/3.0  1 1 0
        1.0    1 1 1 ];
    elseif strcmp(colors, 'copper3');
      dat=[
        0.0  1 1 1
        0.4  0 0 1
        0.6  1 0 1
        0.8  1 0 0
        1.0  1 1 0 ];
    elseif length(colors)>1
      dat=[];
      for i=1:length(colors)
        v=(i-1)/(length(colors)-1);
        if     colors(i) == 'W'; c=[1 1 1];
        elseif colors(i) == 'R'; c=[1 0 0];
        elseif colors(i) == 'G'; c=[0 1 0];
        elseif colors(i) == 'B'; c=[0 0 1];
        elseif colors(i) == 'C'; c=[0 1 1];
        elseif colors(i) == 'M'; c=[1 0 1];
        elseif colors(i) == 'Y'; c=[1 1 0];
        elseif colors(i) == 'K'; c=[0 0 0];
        else
          fprintf('ERROR: bad colors: %s\n', colors);
          return
        end
        dat=[dat; v c];
      end
    else
      fprintf('ERROR: bad colors: %s\n', colors);
      return
    end

    amax=max(max(amp));
    amin=min(min(amp));
    scale=amax-amin;
    
    if swapxy>0
      h = image(time, freq,...
        rb((amp-amin)/scale, dat));
      %xlim([time(1) time(end)]);
      xlim([min(time) max(time)]);
      ylim([min(freq) max(freq)]);
      %ylim([freq(1) freq(end)]);
    else
      h = image(freq, time,...
        rb((amp-amin)/scale, dat));
      %xlim([freq(1) freq(end)]);
      xlim([min(freq) max(freq)]);
      ylim([min(time) max(time)]);
      %ylim([time(1) time(end)]);
    end

    if cbar
      a1=get(h,'Parent');
%      axis([freq(1) freq(end) time(1) time(end)]);
      pos=get(a1, 'Position');
      pos(1)=pos(1)+pos(3); pos(3)=0.001;
      a2=axes('Position', pos, 'YAxisLocation', 'right', 'YLim', [amin0 amax0]);
      data = linspace(amin0,amax0,length(amp(1,:)));
      axes(a1);
    end

%      xlim([time(1) time(end)]);
%      ylim([freq(1) freq(end)]);
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

