function wz = get_wz(hmin, plotting)
%Gives out axial resonance frequency wz for given Imin at 0.5 bar
%Based on fit made from data
%osc2011/old/results/results_wz_20111207-20111213.dat and
%osc2011/HminStates/AllData fullspectra-fits
%Assuming dependence: wz = A*sqrt(Hmin-Hmin0) 
%and w_(mn) = w_L + w_r(m+1) + w_z(n+1/2)

if nargin<2; plotting=0; end

%data from 20111207 and 20111203
Imin = [150 150 200 250 250 350 350 450 450 1500 1500 1250 1000];
w_z  = [25.9 25.2 32.0 35.4 37.2 45.5 48.3 53.4 55.1 104.6 99.1 92.1 83.4];
 
%data from fits in osc2011/HminStates/AllData
[omegas fwids f_z d2 f_r d3 f_r_corr d4 g_z d5 g_r Hmin_0 ~]=...
    textread('/rota/Analysis/PS/osc2011/HminStates/AllData','%f %f %f %s %f %s %f %s %f %s %f %f %[^\n]', 'commentstyle','shell');

f_z = 0.5*f_z;
P = lsqcurvefit(@fitfun,[5 80],Imin,w_z); %P(1) = A, P(2) = Hmin0



if plotting
    find_figure('wz vs Hmin'); clf; hold on;
    h(1)=plot(Imin,w_z,'ob');
    x = linspace(P(2),5000,1000);
    h(2)=plot(x,fitfun(P,x),'-g','LineWidth',2);
end
    %plot also AllData-fits
ave_ind = [];
for i=2:length(f_z)
    if 60 < Hmin_0(i) && Hmin_0(i) < 100
        x = linspace(Hmin_0(i),5000,1000);
        if plotting 
            h(3)=plot(x,fitfun([f_z(i) Hmin_0(i)],x),'-r'); 
        end
        ave_ind = [ave_ind i];
    end
end

f_z_mean = mean([P(1) f_z(ave_ind)']);
Hmin_0_mean = mean([P(2) Hmin_0(ave_ind)']);
if plotting
    x = linspace(Hmin_0_mean,5000,1000);
    h(4)=plot(x,fitfun([f_z_mean Hmin_0_mean],x),'-k','LineWidth',2);
    uistack(h(1),'up',40);
    uistack(h(2),'up',40);

    legend(h,'Data from 12/2011','Fit to data 12/2011','Fits for levels in /HminStates/AllData', ...
        'Averaged fit','Location','NorthWest');
end

wz = fitfun([f_z_mean Hmin_0_mean],hmin);
end

function F = fitfun(p,Imin)
    F = p(1)*sqrt(Imin-p(2));
end

