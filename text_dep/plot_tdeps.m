temps=[0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
omegas=0:0.1:4;

temps=[0.1];
omegas=[5];

dd=zeros(length(temps), length(omegas),1);
rr=zeros(length(temps), length(omegas),1);
aa=zeros(length(temps), length(omegas),1);
bb=zeros(length(temps), length(omegas),1);

ksi=44.2*10^(-6); %dipolar length
c=-48/65*(ksi)^2; %the constant in front of the laplace operator

for i=1:length(temps)
  for j=1:length(omegas)

    unix(sprintf('sed -r -i ''1s/[0-9\\.]+/%.2f/; 7,8s/[0-9\\.]+/%.2f/'' initials.dat',...
       temps(i), omegas(j)));
    [~, ~] = unix('./texture');

    [r a b ~] = textread('texture.dat', '%f %f %f %f',...
                         'commentstyle', 'shell');

    dbdr=(b(2)-b(1))/(r(2)-r(1));

    U1=sin(pi/180*b).^2;
    U2=(pi/180*dbdr*r).^2;
    [en eigv]   = energies_M05bar(r, U1, c);
    [enp eigvp] = energies_M05bar(r, U2, c);
  end
end

figure; clf; hold on;
plot(r, U1, 'r.-');
plot(r, U2, 'b.-');

for i=1:length(en)
  if en(i) > max(U1); break; end;
  plot([0 0.3], en(i)*[1 1], 'r-');
  plot([0 0.3], enp(i)*[1 1], 'b-');
  plot(r(1:end-1), en(i) + 0.1*eigv(:,i)/eigv(1,1), 'r-');
  plot(r(1:end-1), enp(i) + 0.1*eigvp(:,i)/eigvp(1,1), 'b-');
end


ylim([0 max(U1)]);
printf('%f\n', enp(1)/2 );
printf('%f\n', en(1)/2 );


printf('%f\n', (pi/180*dbdr)/4 );

function nemerov(U)
end
