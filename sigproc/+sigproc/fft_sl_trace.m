function [fre, amp, dis] = fft_sl_trace(F, A, f0, df, typ, flock)
% Trace peak position
% Read comments in peak_moments func!
  if nargin<6; flock=10; end
  if nargin<5; typ='min'; end
  for j=1:length(A(1,:))
    [fre(j), amp(j), dis(j)] = sigproc.peak_moments(F, abs(A(:,j))', f0, df,typ,flock);
    if dis(j)~=0; f0=fre(j); else fre(j)=f0; end
  end
end
