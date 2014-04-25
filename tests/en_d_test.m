function en_d_test()
  figure; hold on;

  for a=0:0.5:pi
    for b=0:0.5:pi
      t=0:0.5:pi;
      for i=1:length(t);
        r = abt2r(a,b,t(i));
        E0(i) = en_d0(r);
        E1(i) = en_d1(r);
        E2(i) = en_d2(a,b,t(i));
      end
      plot(t,E1-E0, 'r');
      plot(t,E2-E0, 'b');
    end
  end
end

%% three equivalent expersions for dipolar energy
function e = en_d0(r)
  e=0;
  for i=1:3; for j=1:3;
    e = e + r(i,i)*r(j,j) + r(i,j)*r(j,i);
  end; end
end
function e = en_d1(r)
  e = (trace(r^2) + trace(r)^2);
end
function e = en_d2(a,b,t)
  e = 1/2 * (1 + 4*cos(t))^2 - 0.5;
end
