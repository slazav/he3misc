function en_d_test()
  figure; hold on;

  for a=1:20:180
    for b=1:20:180
      t=1:10:180;
      o=ones(size(t));
      E0 = arrayfun(@en_d_0, a*o,b*o,t);
      E = arrayfun(@en_d, a*o,b*o,t);
      plot(t,E-E0);
%      plot(t,E);
    end
  end
end