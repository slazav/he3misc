function  sliding_fftc(dstr, xfile,tran,fdran,fcorr)
    [tx, x] = read_osc(dstr, xfile);
    L=length(tx);
    
    l = 2^15;
    %l = 25000;
    step = l/8;
    
    fran = [10 2000];
    
    lincur = 1; % fit current with straight line
    
    if nargin < 3, tran = tx([1 end]); end
    if nargin < 5, fcorr = 0; end
    
    i = 1;
    j = 0; %The index determining the time window.
    
    while i < L-l
        tj = tx(i+step); %The recorded times are the end times of each window.
        if tj >= tran(1) && tj <= tran(2)
            j=j+1;
            fprintf('%2d %%', int32(i/L*100));
            [freq_fft, ampj(:,j)] = my_fft(tx(i:i+l-1), x(i:i+l-1),fran);
            time(j) = tj;
            fprintf('\b\b\b\b');
        end
        i = i + step;
    end
    % get set current
    s = regexp(xfile,'_([0-9.]*)s_','tokens','once');
    ttrig = str2num(s{1});
    h = nmrdata_update(dstr);
    rd = nmrdata_time_range(h,ttrig-10,ttrig+90,1);
    d = nmrdata(h,rd,[1 3]);
    ctim = d(:,1)-ttrig;
    cval = d(:,2);
    find_figure('Iset');
    clf
    plot(ctim,cval,'-');
    if lincur 
        jj = find(ctim >= tran(1) & ctim <= tran(2));
        b = regress(cval(jj),[ones(length(jj),1) ctim(jj)]);
        hold on
        plot(ctim([1 end]),b(1)+ctim([1 end])*b(2),'-r');
    end
    
    f0 = get_f0_A(dstr);
    hmin = get_ecoil_cur(dstr,ttrig);
    hlar = 2640.26;
    if lincur
        cur = b(1) + b(2)*(time-(tx(step)-tx(1))/2);
    else
        cur = interp1(d(:,1)-ttrig,d(:,2),time-(tx(step)-tx(1))/2);
    end
        
    
    amp = ampj;
    freq = freq_fft;
    
    amp0 = min(min(amp));
    amp1 = max(max(amp));
    amps = min(amp,amp1/5);
    
    %Let's plot the results
    find_figure(['FFT of signal ' xfile]);
    clf
    [Time,Freq] = meshgrid(time-tran(1),freq);
    [Cur,~] = meshgrid(cur,freq);
    Freq = f0/hlar*(hlar-Cur) - Freq/1000 + fcorr;
    surf(Time,Freq,amps, 'EdgeColor', 'none');
    view(90,-90); %view from the top in 3d figure
    set(gca,'XLim',[0 diff(tran)]);
    if nargin >= 4, set(gca,'YLim',fdran); end
    colormap(hot);
    shading interp
    grid off
    
    if 0
    find_figure(['Raw signal: ' xfile]);
    plot(tx,x);
    
    
    find_figure('fft slice');
    clf
    js = round(length(time/2.5));
    plot(freq_fft,amp(:,js));
    end
end

function [f, SP] = my_fft(t, x, fran)
    L = length(t);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    X = fft(x,NFFT)/L;
    samplerate = 1/(t(2)-t(1));
    f = samplerate/2*linspace(0,1,NFFT/2+1);
    SP = 2*abs(X(1:NFFT/2+1));
    j = find(fran(1)<= f & f<=fran(2)); %Some high enough frequency to limit the number of points
    f=f(j);
    SP=SP(j);
end

