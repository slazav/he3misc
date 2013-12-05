function [res_pr, res_lf, err] = get_fork(date,t0, fork)
%[res_pr, res_lf] = get_fork(date,t0,fork='rotafork1')
%
% get fork width (mHz) a time t0
% res_pr is a value in the previouse point
% res_lf is a result of linear fit between two adjecent points

  if nargin < 3; fork='rotafork1'; end
  fmt='yyyymmdd HH:MM:SS';
  tc = datenum(mktstr(date,t0), fmt);
  t1 = datestr(addtodate(tc, -60, 'minute'), fmt);
  t2 = datestr(addtodate(tc, 30, 'minute'), fmt);

%fprintf('## %s %s\n', t1,t2);
    cmd = ['/rota/programs/bin/getmdata ' fork '_width '...
       '-f "' t1 '" ' '-t "' t2 '" '...
       '-0 "' date ' 00:00:00" -s -q'];

    [s,r] = unix(cmd);
    values=str2num(r);

%fprintf('## %s\n', r);

    res_lf=0;
    res_pr=0;
    err=0;
    if length(values(:,:))==0; return; end
    N=length(values(:,1));
    for i=1:N
      res_pr=values(i,2);
      if i<N
          if (values(i,1) < t0) && (values(i+1,1) > t0)
            res_lf = values(i,2) + (values(i+1,2)-values(i,2))/...
              (values(i+1,1)-values(i,1)) * (t0 - values(i,1));
            err = min(abs(res_lf-values(i,2)), abs(res_lf-values(i+1,2)));
            break
          end
      end
    end
    res_lf=res_lf*1000;
    res_pr=res_pr*1000;
    err = err*1000;
end

function r = mktstr(datestr,t)
    h = floor(t/3600);
    m = floor((t - h*3600)/60);
    s = t - h*3600 - m*60;
    r = sprintf('%s %02d:%02d:%06.3f',datestr,h,m,s);
end
