function ttc = fwid_to_ttc_0bar_fast(fw,fnum,forkpar)

vac_w=0.017; %17mHz vacuum width for fork1, which was used in all 0bar measurements
fw=fw-vac_w;
if nargin < 2
    fnum = 2   ;
elseif isstr(fnum)
    fnum = str2num(fnum(end));
end

if fnum==1
    %a=9359.63;
    a=11700; %new calibration, updated 30.4.2010
else
    
    %a=13134.74
    if nargin<3
    %a=17543; %29 bar, new calibration, updated 30.4.2010
    %a=11460; %0.5 bar, assumed a=pf^4*g, where g is a constant geometrical factor measured@29bar
    a=11200; % 0 bar, a(0)=a(29)*[pf(0)/pf(29)]^4.
    else
    a=forkpar
    end
end
    
ttc = zeros(size(fw));

for i = 1:length(fw)
  if isnan(fw(i))
    ttc(i) = nan;
  else
    ttc(i)=1.7712/log(a/fw(i)); %weak coupling+ gap
  end
end
