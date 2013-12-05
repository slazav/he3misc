function show_fork(zf,fn)
    find_figure(['run ' fn]);
    clf
    pinifit = isfield(zf.ini,'refit');
    if pinifit
        ki = ~cellfun(@isempty,zf.ini.refit);
    else
        ki = zeros(size(zf.ini.ttoday));
    end
    pfinfit = isfield(zf.fin,'refit');
    if pfinfit
        kf = ~cellfun(@isempty,zf.fin.refit);
    else
        kf = zeros(size(zf.fin.ttoday));
    end
    freq0 = floor(min([zf.ini.freq; zf.tdep.freq; zf.fin.freq])/100)*100;
    for xx = {'ini','tdep','fin'}
        zf.(xx{1}).freq = zf.(xx{1}).freq - freq0;
    end
    function doplot(n,fld,ylab)
        subplot(4,1,n);
        plot(zf.ini.ttoday(~ki),zf.ini.(fld)(~ki),'ob');
        hold on
        if pinifit
            plot(zf.ini.ttoday(ki),zf.ini.(fld)(ki),'ob',...
                 'MarkerFaceColor','b',...
                 'ButtonDownFcn',{@show_fit_cb,zf.ini.tlab(ki),...
                                zf.ini.specd(ki),zf.ini.refit(ki)});
        end
        plot(zf.tdep.time,zf.tdep.(fld),'-','Color',[0 0.7 0]);
        plot(zf.fin.ttoday(~kf),zf.fin.(fld)(~kf),'ob');
        if pfinfit
            plot(zf.fin.ttoday(kf),zf.fin.(fld)(kf),'ob',...
                 'MarkerFaceColor','b',...
                 'ButtonDownFcn',{@show_fit_cb,zf.fin.tlab(kf),...
                                zf.fin.specd(kf),zf.fin.refit(kf)});
        end
        fixaxes;
        ylabel(ylab);
    end
    doplot(1,'freq',sprintf('Frequency - %d, Hz',freq0));
    a1 = axis;
    title(fn);
    function fixxax
        a = axis; a(1:2) = a1(1:2); axis(a);
    end
    doplot(2,'width','Width, Hz');
    fixxax;
    doplot(3,'ampl','Amplitude, nA');
    fixxax;
    subplot 414
    plot(zf.tdep.time,zf.tdep.omega,'-r');
    fixaxes; fixxax;
    ylabel('\Omega, rad/s');
    xlabel('time, s');
end

function show_fit_cb(h,eventdata,tlab,specd,fres)
    ax = get(h,'Parent');
    fig = get(ax,'Parent');
    cp = get(ax,'CurrentPoint');
    x0 = cp(1,1);
    y0 = cp(1,2);
    btinf = get(fig,'SelectionType');
    xarr = get(h,'XData');
    yarr = get(h,'YData');
    axx = get(ax,'XLim');
    axy = get(ax,'Ylim');
    dx = (x0-xarr)/(axx(2)-axx(1));
    dy = (y0-yarr)/(axy(2)-axy(1));
    [dmin,j] = min(dx.^2+dy.^2);

    if strcmp(btinf,'extend')
        prevfits = get(fig,'UserData');
        if ~iscell(prevfits), prevfits = {}; end
    else
        prevfits = {};
    end
    prevfits(end+1:end+3) = {specd{j},fres{j},tlab{j}};
    set(fig,'UserData',prevfits);
    showfit('fork fit',prevfits{:});
end

