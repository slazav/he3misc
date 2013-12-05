function xx = sig_reconstruct_ph(tx, time, a0, p2f, p2a, p2p)
% xx(n) =
%   sig_reconstruct_ph(tx(n), time(t), a0(t), p2f(2,t), p2a(2,t), p2p(2,t))
%
% Reconstruct signal in tx points from the result of fit2sin_sl()
%
% -- slazav, feb 2012
  np=length(time); % number of input points
  p1=1; % previous and next point
  p2=1;
  ph=[0,0];
  fprintf('reconstructing signal from sin2fit output:   0 %%');
  for i=1:length(tx)

    pr1=int32(100*(i-1)/length(tx));
    pr2=int32(100*i/length(tx));
    if pr1~=pr2 fprintf('\b\b\b\b\b%3d %%', pr2); end

    % set p1 so that time(p1) < tx(i) < time(p1+1)
    while (p1+1 < np) && (time(p1+1) <= tx(i)); p1 = p1+1; end

    dt=tx(i)-time(p1);
    % result
    xx(i) = a0(p1) + ...
            p2a(p1,1) * sin(2*pi*p2f(p1,1)*dt + p2p(p1,1)) +...
            p2a(p1,2) * sin(2*pi*p2f(p1,2)*dt + p2p(p1,2));

  end
  fprintf('\b\b\b\b\bok   \n');
  xx=xx';
end
