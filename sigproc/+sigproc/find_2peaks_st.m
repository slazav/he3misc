function [p2f, p2a] = find_2peaks_st(tx,xx,window,smooth_win, fmin, fmax)
% [p2f, p2a] = find_2peaks_st(tx,xx,window,smooth_win)
%
% Do fft of the first window of the signal
% and find initial positions for 2 peaks
%
% Do not use smooth_win larger then peak width. 5..7 is good
% for narrow peaks
%
% -- slazav, feb 2012.


  if nargin<5; fmin=10; fmax=10000; end
  [f, a] = sigproc.fft(tx(1:window), xx(1:window),fmin,fmax);
  a=abs(a);
length(a)
  % remove narrow peaks
%  for i=2:length(a)-1;
%    if a(i) > a(i-1) && a(i) > a(i+1)
%      a(i) = max(a(i-1), a(i+1));
%    end
%  end

%find_figure('tmp'); clf; hold on;
%plot(f,a,'r-');

  % smooth fft
  if nargin > 3; a=smooth(a,smooth_win); end

%plot(f,a,'b-');
  [pf, pa] = rel2f.find_peaks(f,a);
%plot(pf,pa,'b.');
  [p2f, p2a] = rel2f.select_2peaks(pf,pa);
%plot(p2f,p2a,'m*');
end
