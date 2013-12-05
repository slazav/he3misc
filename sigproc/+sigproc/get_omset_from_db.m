function [omset,tset] = get_omset_from_db(dstr,tim)

% [omset,tset] = get_omset_from_db(dstr,tim)
% get last set omega value before time
% dstr - date string
% tim - seconds from midnight

db = rota_db;
tq = ['''' dstr '''::date + ' num2str(tim) '* ''1s''::interval'];
query = ['select extract(epoch from (tstamp - ''' dstr ...
	 '''::date)), curstate from xlets.commands_sweep_omega_ae '...
	 'where tstamp < ' tq ' order by tstamp desc limit 1'];

curs = exec(db,query);
curs = fetch(curs);
close(curs);

if length(curs.Data) == 2
    omset = str2num(curs.Data{1,2}); 
    tset = curs.Data{1,1};
else
    omset = []; 
    tset = [];
end
