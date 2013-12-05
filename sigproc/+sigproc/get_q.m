function [A, Q, f0] = get_q(date, AB)
% Fits the Q factor for a given date e.g. '20100526'
% AB parameter is channel name, 'A' or 'B'
% Returns [A, Q, f0]
% A : amplitude at f0
% Q : Q factor
% f0 : resonant frequency

    if (nargin < 2) AB='A'; end

    year = date(1:4);
    month = date(5:6);
    day = date(7:8);
    folder = ['/rota/data/' year '/' month '/' day '/'];
    data = importdata([folder '' date  '-Q' AB '.dat']);
    %firstline = data.textdata{1};
    %drive = str2double(firstline(5:10));
    freq = data.data(:,1)*1000;
    G = data.data(:,2);%/drive;

    modelFun =  @(p,x) ( p(1)./( p(2).*hypot(1-(x.^2./p(3).^2), x./(p(3).*p(2))  )));

    %% calculate moments 
    [A,i]=max(G);
    f0=freq(i);
    df = sqrt(sum((freq-f0).^2 .* G) / sum(G)) / 2;
    Q = f0/df;
    para = [A Q f0];
    para = nlinfit(freq, G, modelFun, para);

    %find_figure('QA'); clf; hold on;
    %plot(freq, G, '.',freq, modelFun(para, freq));

    A=para(1);
    Q=para(2);
    f0=para(3);
end

