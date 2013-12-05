function [p2f, p2a] = select_2peaks(p_f, p_a)
% [p2f(2), p2a(2)] = select_2peaks(p_f(n), p_a(n))
%
% Select two main peaks (with distance >50Hz).
% Use output of find_peaks() as input.
%
% -- slazav, feb 2012.

    [~, max1ind] = max(p_a);

    ii=find(abs(p_f-p_f(max1ind))>10);
    [~, max2ind] = max(p_a(ii));
    max2ind=ii(max2ind);

    p2f = [p_f(max1ind), p_f(max2ind)];
    p2a = [p_a(max1ind), p_a(max2ind)];
    if p2f(1) < p2f(2)
      p2f = [p2f(2), p2f(1)];
      p2a = [p2a(2), p2a(1)];
    end
end
