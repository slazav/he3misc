% Plot spectrum and add  callbacks for peak searching
% Controls (work on both panels):
%   BTN: scroll to new position
%   CTRL-BTN: add point
%   SHIFT-BTN: remove point
%   0-9,+,-: change N_r
%   s - save data to <name>_cwpeaks.mat
%   t - save peak positions to <name>_cwpeaks.txt
%   h - print this help message
% (if <name>_cwpeaks.mat exists then data is readed from it in the beginning)

function plot_spectrum(I, A, name)

  find_figure(name); clf;
  set(gcf, 'KeyPressFcn',@key_cb);

  subplot(2,1,1); hold on;
  plot(I, real(A), 'b-', 'ButtonDownFcn',@btn_cb);
  plot(I, imag(A), 'r-', 'ButtonDownFcn',@btn_cb);
  set(gca, 'ButtonDownFcn',@btn_cba);

  subplot(2,1,2); hold on;
  plot(I, real(A), 'b-', 'ButtonDownFcn',@btn_cb);
  plot(I, imag(A), 'r-', 'ButtonDownFcn',@btn_cb);
  set(gca, 'ButtonDownFcn',@btn_cba);

  title('(press h for help)');

  file = [name '_cwpeaks.mat'];
  if ~unix(['test -f ' file]);
    load(file, 'data');
    data.name = name;
  else
    data.name = name;
    data.NR = [];
    data.I = [];
    data.A = [];
    data.curr_nr = 0;
    data.keep = [2 3]; % keep N objects (2 plots and title...)
    data.xlim = get(gca, 'Xlim');
    data.ylim = get(gca, 'Ylim');
  end
  set(gcf, 'UserData', data);

  set_large_view(mean(data.xlim), 0.5);
  update_marks(1);
end



%% keyboard callback
%%  0..9 set Nz
%%  + inc Nz
%%  - dec Nz
function key_cb(h,eventdata)
  c = eventdata.Character;
  data = get(gcf, 'UserData');

  %% h - print help message
  if c == 'h';
    fprintf('Controls (work on both panels):\n');
    fprintf('  BTN: scroll to new position\n');
    fprintf('  CTRL-BTN: add point\n');
    fprintf('  SHIFT-BTN: remove point\n');
    fprintf('  0-9,+,-: change N_r\n');
    fprintf('  p - update Nr-Nz plot\n');
    fprintf('  s - save data to <name>_cwpeaks.mat\n');
    fprintf('  t - save peak positions to <name>_cwpeaks.txt\n');
    fprintf('  d - save NMR data to <name>_cwdata.txt\n');
    fprintf('  h - print this help message\n');
  end

  %% p - plot Nr-Nz
  if c == 'p';
    h=gcf;
    find_figure('nr, nz'); clf; hold on;
    plot_nr_nz(data.NR, data.I);
    figure(h); % back to main figure
  end

  %% s - save to .mat file
  if c == 's';
    file = [data.name '_cwpeaks.mat'];
    subplot(2,1,2);
    data.xlim=get(gca, 'XLim');
    data.ylim=get(gca, 'YLim');
    save(file, 'data');
    fprintf('Data saved to %s\n', file);
    return;
  end

  %% t - write to text file
  if c == 't';
    file = [data.name '_cwpeaks.txt'];
    fo=fopen(file, 'w');
    fprintf(fo, '# NR NZ Iset\n');
    for i=min(data.NR):max(data.NR);
      ii=find(data.NR==i);
      [~, is] = sort(data.I(ii), 'descend');
      nr=ones(size(ii))*i;
      nz=0:length(ii)-1;
      fprintf(fo, '%02d %02d %f\n', [nr; nz; data.I(is)]);
      fprintf(fo, '\n');
    end
    fclose(fo);
    fprintf('Text data saved to %s\n', file);
    return;
  end

  %% d - write to text file
  if c == 'd';
    file = [data.name '_cwdata.txt'];
    fo=fopen(file, 'w');
    L = get(gca, 'children');
    X = get(L(end-1), 'xdata');
    Y1 = get(L(end-1), 'ydata');
    Y2 = get(L(end), 'ydata');

    fprintf(fo, '# I Abs Disp\n');
    fprintf(fo, '%10f %10e %10e\n', [X; Y1; Y2]);
    fclose(fo);
    fprintf('NMR data saved to %s\n', file);
    return;
  end

  %% +,-,0-9 -- set nr
  nr = str2num(c);
  if c == '+'; nr=data.curr_nr+1; end
  if c == '-'; nr=data.curr_nr-1; end
  if (length(nr)==1);
    data.curr_nr = nr;
    set(gcf, 'UserData', data);
    fprintf('Nr =  %d\n', nr);
    update_marks(0);
  end
end

%% change limits of the lower plot
function set_large_view(x, dI)
  subplot(2,1,2);
%  ylim(get(gca, 'YLim'));
  ylim('auto');
  xlim(x + dI*[-1 1]);
  update_marks(0);
end

%% update marks. if redraw_all == 0 update only view range and Nr= text
function update_marks(redraw_all)
  data = get(gcf, 'UserData');
  col = 'rgbcmk';


 % clear only the last 2 objects
%  subplot(2,1,1);
%  ch = get(gca, 'Children');
%  if length(ch)>=data.keep(1); delete(ch(1:2)); end

%  if (redraw_all)
    for i=1:2
      subplot(2,1,i);

      %% clear all data except plots
      ch = get(gca, 'Children');
      if length(ch)>data.keep(i); delete(ch(1:end-data.keep)); end


      for j=1:length(data.I);
        cn = mod(data.NR(j), length(col))+1;
        plot(data.I(j), data.A(j), ['*' col(cn)], 'ButtonDownFcn',@btn_cb);
      end
    end
%  end

  % plot large view range
  subplot(2,1,2);
  xl =get(gca, 'XLim');
  yl =get(gca, 'YLim');
  subplot(2,1,1);
  plot(xl([1 2 2 1 1]), yl([1 1 2 2 1]), 'k--', 'ButtonDownFcn',@btn_cba);

  txt={['current N_R: ' num2str(data.curr_nr)]};
  text(0.02, 0.1, txt, 'Units','normalized','EdgeColor','red');

  % external plot
  if redraw_all
    h=gcf;
    find_figure('nr, nz'); clf; hold on;
    plot_nr_nz(data.NR, data.I);
    figure(h); % back to main figure
  end

end

%% callback for data plots
function btn_cb(h,eventdata)
  st = get(gcf,'SelectionType');
  cp = get(gca,'CurrentPoint'); I=cp(1,1); A=cp(1,2);
  %fprintf('key pressed: %f %f (%s)\n', I, A, st);

  if strcmp(st,'normal'); % btn, view large
    set_large_view(I, 0.5);
    return
  end

  data = get(gcf, 'UserData');

  if strcmp(st,'extend'); % shift-btn, delete pt
     N=length(data.I);
     [~,i] = min(abs(data.I-I));
     if (length(i)>0)
       ii=find(1:N ~= i); % all points except nearest one
       data.I=data.I(ii);
       data.A=data.A(ii);
       data.NR=data.NR(ii);
     end
  else % ctrl-btn, add point
    data.I(end+1) = I;
    data.A(end+1) = A;
    data.NR(end+1) = data.curr_nr;
  end

  set(gcf, 'UserData', data);

  %% replot marks
  update_marks(1)
end

%% callback for axis: only select zooming point:
function btn_cba(h,eventdata)
  st = get(gcf,'SelectionType');
  cp = get(gca,'CurrentPoint'); I=cp(1,1);
  if strcmp(st,'normal'); % ctrl-btn, view large
    set_large_view(I, 0.5);
  end
end


function nz=get_nz(nr, I)
  nz=zeros(size(nr));
  for i=min(nr):max(nr)
    ii=find(nr==i);
    [~, is] = sort(I(ii), 'descend');
    nz(ii) = is;
  end
end

% plot
function plot_nr_nz(NR, I)

  % fit all axial peaks with 1 or 2 order polynom,
  % save 2 first terms (radial peak position and distance
  % in harmonic approximation)
  peaks=[];
  for i=0:max(NR);
    ii=find(NR==i);
    if ~length(ii); continue; end
    [~, is] = sort(I(ii), 'descend');
    peaks(i+1).nz = 0:length(ii)-1;
    peaks(i+1).iz = I(ii(is));
    order = min(2, length(ii)-1);
    peaks(i+1).pp=polyfit(peaks(i+1).nz, peaks(i+1).iz, order);
    peaks(i+1).nr = i;
    peaks(i+1).ir = peaks(i+1).pp(end);
    if order>0 && i==0;
      peaks(i+1).diz = peaks(i+1).pp(end-1);
    end
  end

  if length(peaks)<1 || length(peaks(1).diz)<1; return; end

  % Fit first 4 radial peaks (found from axial state peaks)
  % with 2-order polynom and find larmor current.
  NR=[peaks.nr];
  IR=[peaks.ir];
  dIZ = peaks(1).diz;

  order = min(2, length(peaks)-1);
  ii=1:min(4, length(peaks)); % fit only 4 first points
  pp1=polyfit(NR(ii), IR(ii), order);
  I00=pp1(end);   % zero peak position from the fit
  % distance between peaks in harmonic approximation:
  if length(pp1)>1; dIR=pp1(end-1); else dIR=0; end
  Ilarm = I00 - dIR/2 - dIZ/4;

  % fit all points
  order = min(5, length(peaks)-1);
  pp2=polyfit(NR, IR, order);

  % plot axial levels
  subplot(1,2,1); hold on; title('axial levels')
  col = 'rgbcmk';
  for i=1:length(peaks);
    NZ = peaks(i).nz;
    IZ = peaks(i).iz;
    NR = i-1;
    IR = peaks(i).ir;
    cn = mod(i-1, length(col))+1; % color name
    plot(NZ, IR-IZ, ['*-' col(cn)]);
  end
  ylabel('I, mA');
  xlabel('N_z');

  % plot radial states
  subplot(1,2,2); hold on; title('radial levels')
  NR=0:length(peaks)-1;
  IR=[peaks.ir];

  for i=1:length(peaks);
    cn = mod(i-1, length(col))+1; % color name
    plot(NR(i), Ilarm-IR(i), ['*' col(cn)]);
  end
  xx=0:0.1:length(IR)-1;

  plot(xx, Ilarm - polyval(pp1, xx), 'k-');
  plot(xx, Ilarm - polyval(pp2, xx), 'k--');
  txt={
    [ 'I_r = ' num2str(-dIR/2) ' mA']
    [ 'I_z = ' num2str(-dIZ/2) ' mA']
    [ 'I_{00} = ' num2str(I00) ]
    [ 'I_{larm} = ' num2str(Ilarm) ]
  };
  text(0.02, 0.90, txt, 'Units','normalized','EdgeColor','red');

  fprintf('%.3f %.3f %.5f %.5f\n', Ilarm, I00, -dIR/2, -dIZ/2);

  ylabel('dF_r, Hz');
  xlabel('N_r');

end


function f = i2f(Iset, f0, Ilarm)
   f  =  f0*(1 - (Iset-Ilarm).*Iset/(Ilarm^2));
end
