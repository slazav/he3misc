
poolsize=7;
if matlabpool('size')==0    
    matlabpool('local', num2str(poolsize))
elseif matlabpool('size')<poolsize
    matlabpool close
    matlabpool('local', num2str(poolsize))
end

TTc=0.13:0.005:0.5;
% 
w0=826*10^3*2*pi;
Press=0.5;

% w0=460*10^3*2*pi;
% P=0;



Diffc=zeros(size(TTc));
parfor k=1:length(TTc)    
   Diffc(k)=diff_coeff(Press,TTc(k),w0); 
end
%%
addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Diffusion coeff/')
PRL65_vs_Lancaster_formula_and_measured_D

plot(h2,TTc,Diffc)

semilogy(h1,1./TTc,Diffc)