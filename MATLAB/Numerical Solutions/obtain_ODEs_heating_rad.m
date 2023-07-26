% ==============================================================================
% This is a function creating the system of ODEs based on the finite difference
% method (FDM) for the ODE solver, during the heating period.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = obtain_ODEs_heating_rad(t,y,A,input)

% Extract some important input data
N = input.N;
dx = input.dx;
Bib = input.Bib;
Rrad = input.Rrad;
SB = input.SB;
lambda = input.lambda1;
Tw = input.Tw;
Tshelf = input.Tshelf(input.time_dim(t));
Tshelf_d = input.temp_non(Tshelf);
phi = input.phi;
temp_dim = input.temp_dim;

% Create a RHS vector
RHS = lambda*ones(N,1);
for i = 1:N
    Qrad = SB*(temp_dim(y(i))^4-Tw^4)/Rrad;
    RHS(i) = RHS(i) - phi(Qrad);
end
RHS(1) = RHS(1);
RHS(N) = RHS(N) + 2*dx*Bib*Tshelf_d/dx^2;

% Create a system of ODEs and export
dydt = A*y(1:end-1) + RHS;
outputs = [dydt; y(end)-Tshelf_d];

return