% ==============================================================================
% This is a top-level routine for view factor calculation.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Numerical Solutions', 'Events', 'View Factor', genpath('View Factor'), genpath('SavedData'));

input_def = get_input_data;  % default inputs
input_def = input_processing(input_def);


%% For obtaining the matrix of view factors
% Change any inputs you want (below listed the frequent ones)
input = input_def;
input.d = 0.01;  % vial diameter (m)
input.c = 0.005;  % distance between vials (m)
input.Lc = 0.3;  % chamber side (m)
input.nx_vial = 10;  % number of vials in x-direction
input.ny_vial = 10;  % number of vials in y-direction
input.layout = [];  % leave empty [] for rectangular or 'CP' for ClosestPacking
input.calmode = 'General';
input = input_processing(input);

figure;
create_VialLayout(input);

vf = get_ViewFactor_network(input);
F = vf.F;

