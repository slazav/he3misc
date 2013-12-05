function [fid, time, omega, fsetA, fmeasA, absA, dispA, fsetB, fmeasB, absB, dispB] = get_nmr(dstr, t1, t2)
% get nmr parameters from signal filename
% example: sigproc.nmr_get('20120705', 48254.32760, 48255.32760)

    fname=['/rota/data/' dstr(1:4) '/' dstr(5:6) '/' dstr(7:8) '/' dstr '.bindat'];
    [st, res] = unix(['/rota/programs/bin/get_nmrtab ' num2str(t1) ' ' num2str(t2) ' ' fname ' | tail -n+2' ]);


    [a, n, e, ~] = sscanf(res, '%f %f %f  %f %f %f %f  %f %f %f %f');
    for i=0:(n-1)/11
      fid(i+1)   = a(i*11+1);
      time(i+1)  = a(i*11+2);
      omega(i+1) = a(i*11+3);

      fsetA(i+1)  = a(i*11+4);
      fmeasA(i+1) = a(i*11+5);
      absA(i+1)   = a(i*11+6);
      dispA(i+1)  = a(i*11+7);

      fsetB(i+1)  = a(i*11+8);
      fmeasB(i+1) = a(i*11+9);
      absB(i+1)   = a(i*11+10);
      dispB(i+1)  = a(i*11+11);
    end
end
