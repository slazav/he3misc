function res = read_data(file, pars)

  % get parameters 
  pp.read_par  = sigproc.par_get('read_par',   pars, 0 ); % read one parameter
  pp.read_freq  = sigproc.par_get('read_freq',   pars, 0 ); % read freq from filename
  pp.read_fpars = sigproc.par_get('read_fpars',   pars, 0 ); % read other parameters from filename
  pp.read_ttrig = sigproc.par_get('read_ttrig',   pars, 0 );
  pp.fork      = sigproc.par_get('fork',       pars, 0 ); % 0, 1, 2
  pp.fork_corr = sigproc.par_get('fork_corr',  pars, 0 ); % fork correction
  pp.max_fork_err = sigproc.par_get('max_fork_err',  pars, 0 );
  pp.max_rel_err  = sigproc.par_get('max_rel_err',   pars, 0 );
  pp.realfre      = sigproc.par_get('realfre',   pars, 0 ); % do not modify freq
  

  if pp.read_par
    [dstr,xfile,a,aerr,a2,t,t2,te,b,be,resnorm,fend,df,tf, par, par1] = ...
      textread([file '.fit1'],...
        '%s %s  %f %f %f %f %f %f %f %f %f %f %f %f  %f %[^\n]','commentstyle','shell'); 
  else
      file
    [dstr,xfile,a,aerr,a2,t,t2,te,b,be,resnorm,fend,df,tf,par1] = ...
      textread([file '.fit1'],...
        '%s %s  %f %f %f %f %f %f %f %f %f %f %f %f %[^\n]','commentstyle','shell'); 
    par=zeros(size(t));
  end
  
  for k=1:length(xfile)
      %str=regexp(par1{k,1},'\s+','split')
      str=par1{k,1};
      ind1=findstr(str,'var=');
      ind2=findstr(str(ind1:end),' ');
      
      if isempty(ind1)
          s=load([file '.cache/fit/' dstr{k} '_' xfile{k} '.mat']);
      else
          aver=str(ind1+4:ind2(1)+ind1-2);
          s=load([file '.cache/fit/' dstr{k} '_' xfile{k} aver '.mat']);
      end
      
      res.residual{k}=s.res.residual;      
      res.timec{k}=s.res.time;
  end

  %% read frequency from filename
%  if pp.read_freq
%    for i=1:length(xfile)
%
%      s = regexp(xfile{i}, '_li([0-9\.]*)Hz','tokens','once');
%      if length(s) > 0;
%        fl(i) = str2num(s{1});
%      else
%        s = regexp(xfile{i}, '_([0-9\.]*)Hz','tokens','once'); % old
%        if length(s) > 0;
%          fl(i) = str2num(s{1});
%        else
%          fl(i)=0;
%        end
%      end
%
%      s = regexp(xfile{i}, '_pi([0-9\.]*)Hz','tokens','once');
%      if length(s) > 0;
%        fp(i) = str2num(s{1});
%      else
%        s = regexp(xfile{i}, '_([0-9\.]*)Hz','tokens','once'); % old
%        if length(s) > 0;
%          fp(i) = str2num(s{1});
%        else
%          fp(i)=0;
%        end
%      end
%    end
%  end
%

%  %% read trigger time from filename
%  if pp.read_ttrig
%    for i=1:length(xfile)
%      s = regexp(xfile{i}, '^([0-9.]*)s','tokens','once');
%
%      if length(s) < 1; s{i} = '0'; end;
%      ttrig(i) = str2num(s{1});
%    end
%    res.ttrig = ttrig';
%  end

  %% read fork if needed
  if pp.fork
    [~,~,fw1p,fw1,fw2p,fw2,~,f1e,f2e, ~] = ...
      textread([file '.frk'],...
        '%s %s  %f %f %f %f  %f  %f %f %[^\n]','commentstyle','shell');

    % check number of signals in the file
    if length(t) ~= length(fw1p)
      fprintf('ERROR: different number of signals in *.fit1 and *.frk\n');
    end

    % select fork 1 or 2
    if pp.fork == 1;
      fw=fw1; fwp=fw1p; fe=f1e;
    else
      fw=fw2; fwp=fw2p; fe=f2e;
    end

    if (fw(end)==0); %if last value not ready
       fprintf('WARNING: last fork value is not ready yet\n');
       fw(end) = fwp(end);
    end
    ii=find(fw==0); fw(ii)=fwp(ii);

    % add fork correction
    fw  = fw - pp.fork_corr;
  else
    fw=zeros(size(t));
    fe=zeros(size(t));
  end

  r=1./t; % relaxation
  r2=1./t2;
  % errors
  te(find(isnan(te))) = 0;
  fe(find(isnan(fe))) = 0;
  re = te./t .* r;

  if pp.max_fork_err == 0; pp.max_fork_err=max(fe./fw)+1; end
  if pp.max_rel_err == 0;  pp.max_rel_err=max(re./r)+1; end


  ii=find(fw==0 | fe./fw <= pp.max_fork_err &...
          re./r  <= pp.max_rel_err);

  res.rel = r(ii);
  res.rel2=r2(ii);
  res.a   = a(ii);
  res.a2=a2(ii);
  res.b   = b(ii);
  res.tf   = tf(ii);
  res.fend  = fend(ii);
  res.rerr  = re(ii);
  res.aerr  = aerr(ii);
  res.resnorm=resnorm(ii);


  if pp.fork
    res.fw  = fw(ii);
    res.fwerr  = fe(ii);
  end
  if pp.read_par
    res.par   = par(ii);
  end

  if pp.read_fpars
    for i = 1:length(ii)
      fpars = sigproc.osc_pars1(xfile{ii(i)});

      s = regexp(par1{ii(i)}, 't0=([0-9.]*)','tokens','once');
      if length(s) < 1; s{1} = '0'; end;
      fpars.time = fpars.time + str2num(s{1});

      res.time(i,1)  = fpars.time;
      res.iset(i,1)  = fpars.iset;
      res.hmin(i,1)  = fpars.hmin;
      res.sens(i,1)  = fpars.sens;
      res.omega(i,1) = fpars.omega;
      res.fl(i,1)    = fpars.fl;
      res.fp(i,1)    = fpars.fp;
      res.dem(i,1)   = fpars.dem;
      if (pp.realfre)
        res.fre(i,1) = res.fend(i);
      else
        res.fre(i,1) = 1000*round(res.fl(i)/1000) - res.fend(i);
      end
    end
  end

end

