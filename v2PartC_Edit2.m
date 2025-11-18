% LMS Part C Code. Poisson-Schrodinger Solver
% Iterative solver for semiconductor quantum well structures

% Preset parameters (Do Not Change)

hbar=1.06e-34; % Plancks constant [Js]
q=1.6e-19;     % Elementary charge [C]
eps0=8.85e-12; % Permitivity of free space [F/m]
epsr=4;        % Relative permittivity
m=0.25*9.1e-31;% Effective mass [kg]
k=8.617e-5;    % Boltzmann constant [eV/K]

% Calculation parameters (Figures need to be inputted)

a=3e-11;       % mesh size [m] (30 pm, DO NOT CHANGE)
mu=   1;       % Chemical potential [eV]
Nout=   round(2e-9/a); % Outer passivating layer (2 nm)
Nw1=    round(3e-9/a); % First (left) well (3 nm) -> 100
Nw2=    round(5.1e-9/a); % Second (right) well (5.1 nm) -> 170
Nin_width_m = 6e-9; % Placeholder: 6 nm
Nin=    round(Nin_width_m/a); % Inner barrier width in mesh units
E_out= 3;      % Outer barrier height
E_m= E_out;    % Inner barrier height set equal to outer barrier height
T=    300;     % Temperature [K]
Vg=0;          % Gate potential (not used for this calculation)

% Preset parameters (Do Not Change)

t0=(hbar^2)/(2*m*(a^2)*q); % Scaling factor for kinetic energy
e0=q*a/eps0;               % Scaling factor for Poisson equation
kT=k*T;                    % Thermal energy in eV
n0=m*kT*q/(2*pi*(hbar^2)); % 2D DOS factor
Np=2*Nout+Nin+Nw1+Nw2; % Total layer thickness (number of mesh points)

% Matrices needed for setting up calculation.

XX=a*1e9*[1:1:Np]; % Position in nm
Ec=[E_out*ones(Nout,1);
    0*ones(Nw1,1);
    E_m*ones(Nin,1);
    0*ones(Nw2,1);
    E_out*ones(Nout,1)];
T=(2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
D2=epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-diag(ones(1,Np-1),-1));
iD2=inv(D2);

Ubdy=-4*[Vg;zeros(Np-2,1);Vg];
U0=iD2*Ubdy;
U1=1e-9*ones(Np,1);
UU=U1;
change=1;
while change>1e-6
    U1=U1+0.1*(UU-U1);
    [P,D]=eig(T+diag(Ec)+diag(U1));
    D=diag(D); % Energy Eigenvalues
    rho=log(1+exp((mu-D)./kT));
    rho=P*diag(rho)*P';
    n=2*n0*diag(rho); 
    UU=U0+iD2*e0*n;
    change=max(max((abs(UU-U1))));
    disp(change)
end;
U=Ec+U1; % Final self-consistent conduction band profile (potential)
ns=1e-4*sum(sum(n.*[zeros(Nout,1);ones(Nw1,1);zeros(Nin,1);ones(Nw2,1);zeros(Nout,1)]));
nn=1e-6*n./a;
disp(['Total sheet density (ns) in cm^-2: ', num2str(ns)])
disp(['Final change (max potential difference): ', num2str(change)])