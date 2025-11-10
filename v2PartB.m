% LMS Part B Code. Calculate energy levels for Finite Well.

% Preset parameters (Do Not Change)

hbar=1.06e-34; %Plancks constant [Js]
q=1.6e-19; %Elementary charge [C]
eps0=8.85e-12; %Permitivity of free space [F/m]
epsr=4; %Relative permittivity
m=0.25*9.1e-31; %Effective mass [kg]
k=8.617e-5; %Boltzmann constant [eV/K]

% Calculation parameters (Change as necessary)
% Default values are:
% a=3e-11m ; mu=1 ; T=300 ; Nout=2e-9/a (Mesh Widths) ; Nw=13nm ; E_out=3eV

a=3e-11;            % mesh size [m] (Default is 30 pm)
mu= 1;              % Fermi level [eV]
T=   300 ;          % Temperature [K]
Nout=  round(2e-9/a)  % Outer passivating layer in mesh units [mesh widths]
Nw=   round(13e-9/a)  % Well width [mesh widths] [nm] 
E_out= 3 ;          % Energy Height of outer barrier (eV)
Vg=0;               % Gate potential (not used for this calculation)

% More parameters needed for calculations. Do not change.

t0=(hbar^2)/(2*m*(a^2)*q); %Scaling factor
e0=q*a/eps0; %Scaling factor
kT=k*T; 
n0=m*kT*q/(2*pi*(hbar^2)); %2D DOS
Np=2*Nout+Nw; %layer thickness in units of mesh size

XX=a*1e9*[1:1:Np];
Ec=[E_out*ones(Nout,1);0*ones(Nw,1);E_out*ones(Nout,1)];
T=(2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
D2=epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-diag(ones(1,Np-1),-1));
iD2=inv(D2);
Vg=0;
Ubdy=-4*[Vg;zeros(Np-2,1);Vg];
U0=iD2*Ubdy;
U1=1e-9*ones(Np,1);UU=U1;change=1;
while change>1e-6
    U1=U1+0.1*(UU-U1);
    [P,D]=eig(T+diag(Ec)+diag(U1));D=diag(D);
    rho=log(1+exp((mu-D)./kT)); rho=P*diag(rho)*P';
    n=2*n0*diag(rho);
    UU=U0+iD2*e0*n;
    change=max(max((abs(UU-U1))));
    U=Ec+U1;
end;
ns=1e-4*sum(sum(n.*[zeros(Nout,1);ones(Nw,1);zeros(Nout,1)]));
nn=1e-6*n./a;
disp('Ground state is')
disp(D(1))