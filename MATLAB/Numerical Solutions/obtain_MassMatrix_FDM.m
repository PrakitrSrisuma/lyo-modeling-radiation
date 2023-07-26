% ==============================================================================
% This is a function creating the mass matrix based on the finite difference
% method (FDM) for the time derivatives in the ODEs.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = obtain_MassMatrix_FDM(input)

% Extract some important input data
N = input.N;

% Create the mass matrix 
M = eye(N);

% Add a dimension for collecting the shelf temperature
M = [M, zeros(N,1)];
M = [M; zeros(1,N+1)];

% Export
outputs = M;

return