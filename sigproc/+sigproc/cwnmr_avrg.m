function [I0 A0 N] = cwnmr_avrg(dstr, t1, t2, pars)
% average all sweeps in t1..t2 interval
%   pars:
% dir    --  select sweep direction, 1 or -1 ( default -1)
% minp   -- set min number of points in selected sweeps (100)
% mint   -- set min time of selected sweeps (100)
% smooth -- set smooth length for base calculation (200)
% chan   -- NMR channel, A or B (A)
%     slazav 07/12/2012

  % get parameters:
  if nargin < 4; pars=''; end
  pp.dir     = sigproc.par_get('dir',    pars, -1 );
  pp.minp    = sigproc.par_get('minp',   pars, 100 );
  pp.mint    = sigproc.par_get('mint',   pars, 100 );
  pp.smooth  = sigproc.par_get('smooth', pars, 200 );
  pp.chan    = sigproc.par_get('chan',   pars, 'A' );

  % get nmr data
  [fid, time, omega, fsetA, fmeasA, absA, dispA, fsetB, fmeasB, absB, dispB] =...
     sigproc.cwnmr_get(dstr, t1, t2);


  pdir=0; % previous sweep direction. 
  ppt=1;  % previous point of sweep direction change. 
  A0=[];  % arrays to collect data
  I0=[];
  N=0;    % number of sweeps

  % select channel
  if pp.chan == 'A';
    fset=fsetA; Abs=absA; Disp=dispA;
  elseif pp.chan == 'BA'
    fset=fsetB; Abs=absA; Disp=dispA;
  elseif pp.chan == 'AB';
    fset=fsetA; Abs=absB; Disp=dispB;
  else
    fset=fsetB; Abs=absB; Disp=dispB;
  end

  % collect all data in I0 A0
  d=0;
  for i=2:length(time)
    % sweep direction in the point.
    if     fset(i) > fset(i-1); d=1;
    elseif fset(i) < fset(i-1); d=-1;
    end

    % process sweep until direction changed or last point reached
    if d == pdir && i ~= length(time); continue; end

    % select sweep according to pars
    if (pdir == pp.dir) &&...                % select by direction
       (time(i)-time(ppt) >= pp.mint) &&...  % select by min time
       (i-ppt+1 >= pp.minp)                  % select by min point number

      A=Abs(ppt:i) + j*Disp(ppt:i);
      I=fset(ppt:i);

      % we need sort+unique for interp1 function
      [I, i1] = sort(I); A=A(i1);
      [I, i1] = unique(I); A=A(i1);

      if length(I0)==0; % fill A0/I0 with the first sweep
        A0=A;
        I0=I;
        N=1;
      else % add another sweeps to A0/I0
        i1  = find(I  >= min(I0) & I  <= max(I0));
        i2  = find(I0 >= min(I)  & I0 <= max(I));
        A0 = A0(i2) + interp1(I(i1),A(i1), I0(i2));
        I0 = I0(i2);
        N=N+1;
      end
    end
    pdir=d;
    ppt=i-1;
  end
  A0 = A0/N;

  if (length(A0)==0) error ('No data in this time range!'); else

  % subtract smoothed base
  A0  = A0 - smooth(real(A0),pp.smooth)' - j*smooth(imag(A0),pp.smooth)';
%  A0  = smooth(real(A0),pp.smooth)' + j*smooth(imag(A0),pp.smooth)';

end
