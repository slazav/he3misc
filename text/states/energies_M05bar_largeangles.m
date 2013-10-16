function [E,Eigenvectors]=energies_M05bar_largeangles(N,T,f_larmor,omega,omegav,lambdaHV,chi,APsi,Initial_texture)
format long;
R_0=0.3*10^(-2);     %Real cylinder radius, meters
R=R_0;  % SIGNIFFICANT cylinder radius, meters. Less than of equal to real radius R_0
ksi=10*10^(-6); %dipolar length
c=-8/15*(ksi)^2; %the constant in front of the laplace operator
k=10/13;

h=R/(N+1);  % h=step size
h_0=R_0/(N+1); %step size in temporary potential
r= (h_0:h:((N-1)*h+h_0))';   %discretization points starting from h_0
r_0= (1:N)'*h_0;%Discretization points in temporary potential



%Import potential here
%NN=N;
%if N>1000-1
%    NN=1000-1;
%end
[textur,~]=textureM05bar(N+1,[T 0.5 f_larmor/1000 0.3 omega omegav -1 lambdaHV chi],500,[0.7 10],Initial_texture,APsi);  %N points = N+1 intervals
%plot the texture
%plot(textur(:,1),textur(:,3),'-')

ttV=(sind(textur(:,3))).^2;

%Interpolating the potential within the signifficant radius
%tV=interp1(r_0,ttV(2:(N+1)),r);   %Length of the potential vector ttv is N+1, ttV(1)=0

V=spdiags(ttV,0,N,N);


%Creating R_L
alpha=textur(:,2);
beta=textur(:,3);
R_Ltemp=zeros(N,1);
for kk=1:N+2
    R_Ltemp(kk)=1/32*(171+35*cosd(2*alpha(kk))+50*(cosd(alpha(kk)))^2*cosd(4*beta(kk))+120*cosd(2*beta(kk))*(sind(alpha(kk)))^2+80*sqrt(15)*cosd(beta(kk))*sind(2*alpha(kk))*(sind(beta(kk)))^2);
end
R_L=spdiags(R_Ltemp(2:N+1),0,N,N);

%Creating the discretization matrix L for the second derivative: diagonal elements are "-2" and on the
%left and on the right side of the twos there are "1":s

tD = sparse(1:N,1:N,-2*ones(1,N),N,N);
tE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
L = tE+tD+tE';


% Boundary conditions for L

    %We demand S'(0)=0 (symmetry reasons), S'(R)=S'(1)=0
    L(1,1)=-2/3;
    L(1,2)=2/3;
    L(N,N)=-2/3;
    L(N,N-1)=2/3;
    
    
L=L*1/(h^2);     %correct scale

%Creating the discretization matrix D for the first derivative : zeros on
%the diagonal, "1" right from the diagonal and "-1" left from it

ttE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
D=-ttE+ttE';
tttE = sparse(2:N+2,1:(N+2)-1,ones(1,(N+2)-1),N+2,N+2);
Dtexture=-tttE+tttE';
%Boundary conditions for D
    
    %We demand S(R)=S(1)=0, again no modification needed
    %We want S'(0)=0, thus
    D(1,1)=-2/3;
    D(1,2)=2/3;
    D(N,N)=-2/3;
    D(N,N-1)=2/3;
    
D=D*1/(2*h); %correct scale
Dtexture=Dtexture*1/(2*h);
%for the radial laplace we need 1/r * d/dr * r * d/dr , not d/dr. Creating a diagonal matrix T of
% "1/r":s and "r":s

T=spdiags(r.^(-1),0,N,N);
%T2=spdiags(r,0,N,N);


%Building the complete discr. matrix K
Identity=spdiags(ones(N,1),0,N,N);
RLdifftemp=Dtexture*R_Ltemp;
R_Ldiff=spdiags(RLdifftemp(2:N+1),0,N,N);
%R_Ldiff(1,1)=R_Ldiff(2,2);
%R_Ldiff(N,N)=R_Ldiff(N-1,N-1);
K=(V+(((1+k)*Identity-k/16*R_L)*(L+T*D)-k/16*R_Ldiff*D).*c); % c: see for the top



%Organizing s (<=N) first eigenvectors into matrix "Eigenvectors" and corresponding eigenvectors in vector "E" 
s=round(N/5);

%Starting vector for eigs
%opts.v0=1/8*exp(-100/R^2*r.^2);
%Eigenvalues
[Eigenvectors,Energies,flag]=eigs(K,s,'SR');%,opts);
while flag ~= 0
    [Eigenvectors,Energies,flag]=eigs(K,s,'SR');
end
%Turning energies into a vector
E=zeros(1,s);
for n=1:s
    E(n)=Energies(n,n);
end
    
end