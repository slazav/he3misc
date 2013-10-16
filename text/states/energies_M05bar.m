function [E,Eigenvectors]=energies_M05bar(N,T,f_larmor,lambda_omega,Omega,Omegav,chi,APsi,Initial_texture)
format long;
R_0=0.3*10^(-2);     %Real cylinder radius, meters
R=R_0;  % DO NOT ADJUST
ksi=44.2*10^(-6);%dipolar length
c=-48/65*(ksi)^2; %the constant in front of the laplace operator


h=R/(N+1);  % h=step size
h_0=R_0/(N+1); %step size in temporary potential
r= (h_0:h:((N-1)*h+h_0))';   %discretization points starting from h_0
r_0= (1:N)'*h_0;%Discretization points in temporary potential







%Import potential here
%NN=N;
%if N>1000-1
%    NN=1000-1;
%end
[textur,~]=textureM05bar(N+1,[T 0.5 f_larmor 0.3 Omega Omegav lambda_omega 0.9 chi],500,[0.7 10],Initial_texture,APsi);  %N points = N+1 intervals
%plot the texture
%plot(textur(:,1),textur(:,3),'-')

ttV=(sind(textur(:,3))).^2;
%Interpolating the potential within the signifficant radius
%tV=interp1(r_0,ttV(2:(N+1)),r);   %Length of the potential vector ttv is N+1, ttV(1)=0

V=spdiags(ttV(2:N+1),0,N,N);

%Creating the discretization matrix L for the second derivative: diagonal elements are "-2" and on the
%left and on the right side of the twos there are "1":s

tD = sparse(1:N,1:N,-2*ones(1,N),N,N);
tE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
L = tE+tD+tE';


% Boundary conditions for L

    %We demand S'(0)=0 (symmetry reasons)
    L(1,1)=-2/3;
    L(1,2)=2/3;
    %On the boundary S(R)=S(1)=0, no modification needed
    
L=L*1/(h^2);     %correct scale



%Creating the discretization matrix D for the first derivative : zeros on
%the diagonal, "1" right from the diagonal and "-1" left from it

ttE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
D=-ttE+ttE';

%Boundary conditions for D
    
    %We demand S(R)=S(1)=0, again no modification needed
    %We want S'(0)=0, thus
    D(1,1)=-2/3;
    D(1,2)=2/3;
   
    
D=D*1/(2*h); %correct scale

%for the radial laplace we need 1/r * d/dr , not d/dr. Creating a diagonal matrix T of
% "1/r" :s

T=spdiags(r.^(-1),0,N,N);


%Building the complete discr. matrix K

K=(V+(L+T*D).*c); % c: see for the top



%Organizing s (<=N) first eigenvectors into matrix "Eigenvectors" and corresponding eigenvectors in vector "E" 
s=round(N/10);


if 1==0 %method 1: fast but has convergence problems sometimes    
    
[Eigenvectors,Energies,flag]=eigs(K,s,'SR');%,opts);
while flag ~= 0
    [Eigenvectors,Energies,flag]=eigs(K,s,'SR');
end
E=zeros(1,s);
for n=1:s
    E(n)=Energies(n,n);
end


else %method 2: no convergence problems, slow
[Eigenvectors,Energies]=eig(full(K));
[E,I]=sort(diag(Energies));
E=E';
Eigenvectors=Eigenvectors(:,I);

end
end