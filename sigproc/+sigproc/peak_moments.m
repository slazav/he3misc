function [fre, amp, dis] = peak_moments(F, A, f0, df, typ, flock)
% Find amplitude, frequency and width of peak near f0 freq.
%
% Amplitude is calculated as a sqrt of first moment of A^2 (power)
% in the f0-df..f0+df range. Initial f0 value is moved to the nearest
% local maximum of smoothed signal
%
% Frequency is calculated depending on typ parameter:
%  'max' (default) - Quadratic fit near maximum. Works better for noisy or
%                    non-symmetric peaks, if non-symmetric side peaks exists.
%                    Works better for very narrow signals (with width approx
%                    equal to fft frequency step).
%  'moment'        - second moment of A^2. Works better for large
%                    lorenzian-style signals.
%
% Dispersion is calculated as third central moment of A^2. Works bad
% for very narrow signals (problem of fft - depandence on phase).
% It is set to 0 if peak is lost (strange behaviour).
%
%    if nargin<6; flock=10; end
%    if nargin<5; typ='max'; end

    P=abs(A).^2;

    %% Lock to line. Go from f0 to the nearest local maximum of smoothed A^2.
    [~, i0] = min(abs(F-f0)); % index of f0:
    i0=0;
    for k=0:100 % go to local maximum
      ii=find(abs(F-f0)<flock);
      [~, j] = max(smooth(P(ii),4)); j=ii(j);
      if i0==j; break; end
      i0=j;
      f0=F(j);
    end

    %% set integrating range
    ii  = find(abs(F-f0)<df);
    if length(ii)<2
      amp=0; fre=0; dis=0; return;
    end

    %% calculate moments
    power=P(ii) .* blackman(length(ii))';

    m0 = sum(power);
    m1 = sum(F(ii) .* power);
    m2 = sum((F(ii)-m1/m0).^2 .* power);

    amp = sqrt(m0);
    fre = m1/m0;
    dis = sqrt(m2/m0);

    if strcmp(typ,'max'); [fre, ~]= find_max(F,P,i0); end

    % "bad" line
    win=8;
    sm0 = sm(P, i0, win);
    if sm0 < (sm(P, ii(1), win) + sm(P, ii(end), win))*0.8;
      dis=0; end

end


function s=sm(arr, i0, win)
  s=0;
  bl=blackman(2*win+1)';
  for i=i0-win:i0+win
    if i > 1 && i <= length(arr)
       s = s + arr(i)*bl(i-i0+win+1);
    end
  end
end

function [f0,a0] = find_max(freq, amp, m)
    if m<2; f0=freq(1); a0=amp(1); return; end
    if m>length(freq)-1; f0=freq(end); a0=amp(end); return; end
    % fit amp(freq) near maximum by Ax^2+Bx+C.
    f1 = freq(m-1);
    f2 = freq(m);
    f3 = freq(m+1);

    a1 = amp(m-1);
    a2 = amp(m);
    a3 = amp(m+1);

    if (a2<a1) || (a2<=a3); 
      f0=freq(m); a0=amp(m); return;
    end

    AA = ((a2-a1)/(f2-f1)-(a3-a1)/(f3-f1))/...
         ((f2^2-f1^2)/(f2-f1)-(f3^2-f1^2)/(f3-f1));
    BB = (a2-a1)/(f2-f1) - AA*(f2^2-f1^2)/(f2-f1);
    CC = a1 - AA * f1^2 - BB * f1;

    f0 = -BB/2.0/AA;
    a0 = AA*f0^2 + BB*f0 + CC;
end
