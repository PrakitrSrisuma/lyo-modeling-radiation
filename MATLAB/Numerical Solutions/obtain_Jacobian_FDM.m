% ==============================================================================
% This is a function creating the Jacobian based on the finite difference
% method (FDM) for the system of ODEs.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = obtain_Jacobian_FDM(input)

% Extract some important input data
N = input.N;
dx = input.dx;
Bib = input.Bib;

% Create a matrix
diagonal = (-2/dx^2)*ones(N,1);
offdiag = (1/dx^2)*ones(N-1,1);
A = diag(diagonal,0) + diag(offdiag,1) + diag(offdiag,-1);

% Stamp some boundary nodes
A(1,1) = -2/dx^2;
A(1,2) = 2/dx^2;
A(N,N) = -(2*dx*Bib+2)/dx^2;
A(N,N-1) = 2/dx^2;

% Export
outputs = A;

return