function en_dr_test()
  figure; hold on;

  for a=0:1:pi
    for b=0:1:pi

      for thx=-0.1:0.1:0.1
      for thy=-0.1:0.1:0.1
      for thz=-0.1:0.1:0.1
      th=[thx thy thz];
      t=0:0.5:pi;

      for i=1:length(t);
        r = abt2r(a,b,t(i));
%        E0(i) = en_dr0(r,th);
        E0(i) = en_dr1(r,th);
        E1(i) = en_dr2(a,b,t(i),th);
      end
      plot(t,E1-E0, 'r');

      end
      end
      end
    end
  end

end

%% dipolar energy
function e = en_d0(r)
  e=0;
  for i=1:3; for j=1:3;
    e = e + r(i,i)*r(j,j) + r(i,j)*r(j,i);
  end; end
end

%% small rotation of matrix -- my
function e = en_dr0(r, th)
  et = [    0   th(3) -th(2)
         -th(3)    0   th(1)
          th(2) -th(1)    0];
  tt = th*th'; %'
  r1 = r;
  for i=1:3; for j=1:3; for k=1:3;
          r1(i,j) = r1(i,j) - et(i,k)*r(k,j) + th(i)*th(k)/2 * r(k,j) - tt/2*r(i,j);
  end;end;end;
  e = en_d0(r1);
end
function e = en_dr1(r, th)
  e=0;
  et = [    0   th(3) -th(2)
         -th(3)    0   th(1)
          th(2) -th(1)    0];
  dd = [1 0 0
        0 1 0
        0 0 1];
  tt = th*th'; %'

  for i=1:3; for j=1:3;
    e = e + (r(i,i)*r(j,j) + r(i,j)*r(j,i))*(1-tt);
    for b=1:3;
      e = e - 2* et(j,b) * r(b,j)*r(i,i);
      e = e - 2* et(j,b) * r(b,i)*r(i,j);

      e = e + th(j)*th(b) * r(b,j)*r(i,i);
      e = e + th(j)*th(b) * r(b,i)*r(i,j);

      for c=1:3;
         e = e + et(i,c)*et(j,b) * r(b,j)*r(c,i);
         e = e + et(i,b)*et(j,c) * r(b,j)*r(c,i);
      end;
    end;
  end; end
%  e = - 4*(4*cos(t)+1) * sin(t) * nt ...
end

function e = en_dr2(a,b,t, th)
  n(1) = sin(b)*cos(a);
  n(2) = sin(b)*sin(a);
  n(3) = cos(b);
  ct=cos(t);
  st=sin(t);
  nt=n*th'; tt=th*th';

  e = - 0.5 + 1/2*(4*ct+1)^2 ...
      - 4*(4*ct+1) * st * nt ...
      - (4*ct^2+5*ct+1)*tt ...
      + (9 + 3*ct - 12*ct^2) *nt^2;
end

