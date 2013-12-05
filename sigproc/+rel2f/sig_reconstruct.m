function xx = sig_reconstruct(tx, time, p2f, p2a)
% xx(n) = sig_reconstruct(tx(n), time(t), p2f(t), p2a(t))
%
% Reconstruct signal for time values xt from the result
% of select_2peaks_t(), frequencies and amplitudes of
% two peaks, p2f(time) and p2a(time)
%
% -- slazav, feb 2012 
  np=length(time); % number of input points
  p1=1; % previous and next point
  p2=1;
  ph=[0,0];
  fprintf('reconstructing signal:   0 %%');
  for i=1:length(tx)
    pr1=int32(100*(i-1)/length(tx));
    pr2=int32(100*i/length(tx));
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); end

    % set p1 and p2 so that time(p1) < tx(i) < time(p2)
    while (p1+1 < np) && (time(p1+1) < tx(i)); p1 = p1+1; end
    p2=p1+1;
    if i>1; dt=tx(i)-tx(i-1); else dt=0; end
    % calculate amplitudes, freqs and phases of two peaks at tx(i)
    for n=1:2
      if time(p2) ~= time(p1)
        a(n) = p2a(n,p1) + (p2a(n,p2) - p2a(n,p1))/(time(p2)-time(p1)) * (tx(i)-time(p1));
        f(n) = p2f(n,p1) + (p2f(n,p2) - p2f(n,p1))/(time(p2)-time(p1)) * (tx(i)-time(p1));
      else
        a(n) = p2a(n,p1);
        f(n) = p2f(n,p1);
      end
      ph(n) = ph(n) + 2*pi*f(n) * dt;
    end
    % result
    xx(i) = a(1) * sin(ph(1)) + a(2) * sin(ph(2));
  end
  fprintf('\b\b\b\b\bok   \n');
  xx=xx';
end
