function [res,ov] = procfork2(fork,dstr,tx)

if nargin == 2
    tx = dstr;
    dstr = fork;
    fork = 'rotafork2';
end

addpath /rota/Analysis/Fork

hdata = nmrdata_update(dstr);

for xx = {'demag_E','demag_K','demag_L','low_T'}
  fdat = load_fork_data([dstr '-' fork '-' xx{1} '.dat']);
  if ~isempty(fdat), break; end
end

ov = [nan nan];
res = [];


% locate and load fork saved data
datad = ['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/'];
dl = dir([datad '/' dstr 'T*-' fork '-*.dat']);
if strcmp(fork,'rotafork2')
    dl = [dl; dir([datad '/' dstr 'T*-rotafork-tdep.dat'])];
end
if isempty(dl), return; end

% calculate times for each file
s = cell2mat(cellfun(@(x)(x(10:15)),{dl.name},'UniformOutput',0));
hms = sscanf(s,'%02d%02d%02d',[3 inf])';
[dtims,ix] = sort(hms(:,1)*3600+hms(:,2)*60+hms(:,3));
dl = dl(ix);
istdep = cellfun(@(x)(strcmp(x(end-7:end),'tdep.dat')),{dl.name})';
jt = find(istdep & dtims <= tx,1,'last');
if isempty(jt), return; end

  fprintf('date: %s\n', dstr);
  fprintf('timedep fork files:');
  fprintf(' %d', dtims);
  fprintf('\n');

% get couple of sweeps in the beginning and at the end
jini = find(~istdep & dtims < tx,2,'last');
if length(jini)==2 && diff(jini)~=1, jini = jini(2); end
t1 = min(dtims([jini;jt])) - 1;

jfin = find(~istdep & dtims > tx,2);
if length(jfin)==2 && diff(jfin)~=1, jfin = jfin(1); end
if ~isempty(jfin) && jfin(end) < length(dtims)
    t4 = dtims(jfin(end)+1);
else
    t4 = 2*90000;
end

  fprintf('time between full sweeps: %d .. %d\n', t1, t4);

if isempty(jini) && isempty(jfin), return; end


fd = load_fork_from_db(fork,dstr,t1,t4);

specini = {}; specfin = {};
tlabini = {}; tlabfin = {};

tdep = [];

t2 = tx;
t3 = t2;

useexc = nan;
% select ini, tdep and fin parts
for j = find(dtims >= t1 & dtims < t4)'
    fn = dl(j).name;
    tlab = fn(1:15);
    if strcmp(fn(end-7:end),'tdep.dat')
        [d,exc,useexc] = readftdep(datad,fn,useexc);
        if isempty(d), continue; end
        if d(1,1) >= t1 && d(end,1) < t4
          % clean possible relaxing part in the beginning
          if size(d,1) >= 2
            if d(2,3) > d(1,3)
              j = d(2:end,3) <= d(1:end-1,3);
            else
              j = d(2:end,3) >= d(1:end-1,3);
            end
            k = find(j,1);
            if ~isempty(k),
              d = d(k:end,:);
            end
          end
            tdep = [tdep; d];
            %fprintf('   %s tdep: %s\n',fork,tlab);
        end
    else
        [d,exc,useexc] = readfswp(datad,fn,useexc);
        if isempty(d), continue; end
        if d(1,1) >= t1 && d(end,1) < t2
            specini{end+1} = d;
            tlabini{end+1} = tlab;
            %fprintf('   %s ini: %s\n',fork,tlab);
        elseif d(1,1) >= t3 && d(end,1) < t4
            specfin{end+1} = d;
            tlabfin{end+1} = tlab;
            %fprintf('   %s fin: %s\n',fork,tlab);
        end
    end
end
% find min and max freq for background type selection
alld = vertcat(specini{:},tdep,specfin{:});
if ~isempty(alld)
    swpwid = max(alld(:,2)) - min(alld(:,2));
    if swpwid < 75
        bgorder = 0; % constant
    elseif swpwid < 500
        bgorder = 1; % linear
    else
        bgorder = 2; % quadratic
    end
    %bgorder
    [cnv,frini,frfin] = convert_tdep(tdep,specini,specfin,bgorder);
else
    cnv = zeros(0,5);
    frini = {};
    frfin = {};
end

res.ini = merge_fdat(cnv(1:length(tlabini),:),...
                     specini,frini,tlabini,fdat,fd,t1,t2);
res.fin = merge_fdat(cnv(end-length(tlabfin)+1:end,:),...
                     specfin,frfin,tlabfin,fdat,fd,t3,t4);
            
cnv = cnv(cnv(:,5)==0,1:4);
% clean relaxing part in the beginning
res.tdep = struct('time',cnv(:,1),'freq',cnv(:,2),...
                  'width',cnv(:,3),'ampl',cnv(:,4),'omega',zeros(0,1));

if ~isempty(cnv)
    omd = nmrdata(hdata,nmrdata_time_range(hdata,cnv(1,1),cnv(end,1),1),[1 2]);
%    res.tdep.omega = interp1(omd(:,1),omd(:,2),cnv(:,1))/1000;
end

if ~isempty(res.ini.width), ov(1) = mean(res.ini.width); end
if ~isempty(res.fin.width), ov(2) = mean(res.fin.width); end

function [d,exc,useexco] = readfswp(datad,fn,useexc)
fi = fopen([datad '/' fn]);
a = fscanf(fi,['# Sweep from %fHz to %fHz step %fHz delay %fs\n'...
               '# Excitation %fVpp attenuator %f dB\n'...
               '# time(s), freq(Hz), absorption(nA), dispersion(nA)\n']);
if length(a) == 6
    exc = a(5);
    [d,useexco] = readfrest(fi,exc,useexc);
else
    exc = nan;
    d = [];
    useexco = useexc;
end
fclose(fi);

function [d,useexco] = readfrest(fi,exc,useexc)
d = fscanf(fi,'%f %f %f %f\n',[4 inf])';
if length(d) ~= 0 && d(end,4) == 0 % partial line
    d(end,:) = [];
end
useexco = useexc;
if ~isempty(d)
    if isnan(useexc)
        useexco = exc;
    elseif exc~=useexc
        d(:,[3 4]) = d(:,[3 4])/exc*useexc;
    end
end

function [d,exc,useexco] = readftdep(datad,fn,useexc)
d = [];
exc = nan;
useexco = useexc;
fi = fopen([datad '/' fn]);
s = fgets(fi);
if length(s) > 4
    if strcmp(s(1:4),'# ex')
        exc = sscanf(s,'# excitation %f');
        nskip = 4;
    else
        exc = sscanf(s,'# Excitation %f');
        nskip = 2;
    end
    for i = 1:nskip
        s = fgets(fi);
    end
    if length(s) > 0
        [d,useexco] = readfrest(fi,exc,useexc);
    end
end
fclose(fi);

function res = merge_fdat(cnv,specd,fres,tlabs,fdat,fd,t1,t2)
% merge fork data from different sources

if ~isempty(fdat)
    % fdat for t1..t2
    j = fdat.ttoday >= t1 & fdat.ttoday <= t2;
    nms = fieldnames(fdat);
    for i = 1:length(nms)
        nm = nms{i};
        fdsel.(nm) = fdat.(nm)(j);
    end
    fdat = fdsel;
end

if isempty(cnv) && isempty(fdat) % data only in database
    j = find(fd(:,1) >= t1 & fd(:,1) <= t2);
    res.ttoday = fd(j,1);
    res.width = fd(j,2);
    res.freq = fd(j,3);
    res.ampl = nan(length(j),1);
    return
end

% fill fdat from the database to overcome precision problem
nfdat = length(fdat.tlab);
for j = 1:nfdat
    k = find(abs(fd(:,1) - fdat.ttoday(j)) < 1, 1);
    if isempty(k), continue; end
    fdat.width(j) = fd(k,2);
    fdat.freq(j) = fd(k,3);
end

% now merge fit results with fdat
res = fdat;
res.refit = cell(nfdat,1);
res.specd = cell(nfdat,1);
for i = 1:length(tlabs)
    tlab = tlabs{i};
    j = find(strcmp(res.tlab,tlab),1);
    if isempty(j)
        j = length(res.tlab)+1;
        res.tlab{j} = tlab;
    end
    res.ttoday(j) = cnv(i,1);
    res.freq(j) = cnv(i,2);
    res.width(j) = cnv(i,3);
    res.ampl(j) = cnv(i,4);
    res.refit{j} = fres{i};
    res.specd{j} = specd{i};
end
