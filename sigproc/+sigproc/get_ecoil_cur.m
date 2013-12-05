function res = get_ecoil_cur(dstr,tim)
% get last set current in the extra coil

  db = rota_db;
  tq = ['''' dstr '''::date + ' num2str(tim) '* ''1s''::interval'];
  %query = sprintf(['select "fsetA" from data.sweep_nmra where tsec ' ...
  %                 '<= from_tstamp(''%s'')+%.3f order by tsec desc limit 1'],dstr,tim);
  query = ['select curstate from xlets.commands_sweep_nmrA where ' ...
            'tstamp <= ' tq ' order by tstamp desc limit 1'];
  curs = exec(db,query);
  curs = fetch(curs);
  res = cell2mat(curs.Data);
  close(curs);
