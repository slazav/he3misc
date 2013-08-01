
poolsize=7;
if matlabpool('size')==0    
    matlabpool('local', num2str(poolsize))
elseif matlabpool('size')<poolsize
    matlabpool close
    matlabpool('local', num2str(poolsize))
end

TTc=0.25;
% 
% w0=826*10^3*2*pi;
% P=0.5;

w0=[400:100:2000]*10^3*2*pi;
Press=0.5;
TTc=0.25*ones(size(w0));


Diffc=zeros(size(TTc));
parfor k=1:length(TTc)    
    
   Diffc(k)=diff_coeff(Press,TTc(k),w0(k)); 
end
%%
figure;
loglog(w0,Diffc);