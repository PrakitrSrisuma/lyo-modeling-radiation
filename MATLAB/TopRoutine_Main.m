% ==============================================================================
% This is a top-level routine.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Numerical Solutions', 'Events', genpath('View Factor'), genpath('Saved Data'));

% Define and extract important input data
input_def = get_input_data;  % default inputs
input_def = input_processing(input_def);

% Option selection
Original = 'off';  % 'on' for normal simulation / no radiation
Sim_Single = 'off';  % 'on' for simplified network for one vial
Sim_Multiple = 'off';  % 'on' simplified network for multiple vials
Network = 'on';  % 'on' for radiation network


%% Default Simulation - No Radiation
switch Original
case 'on'

input = input_def;
input.rad = 'off';
input.mode = 'CFD';
input = input_processing(input);
[t_FDM, t_FDM2, y_FDM, y_FDM2] = obtain_sol_rad(input);
[time_FDM, Temp_FDM, Tshelf_FDM, ~] = temperature_postprocessing(t_FDM, y_FDM, input);
[time_FDM2, itf_FDM, Tshelf_FDM2, tdry_FDM, Tb_FDM] = interface_postprocessingV2(t_FDM2, y_FDM2, input);
tdry_normal = tdry_FDM;

end


%% Independent Radiation
switch Sim_Single
case 'on'

input = input_def;
input.rad = 'sim';
input.mode = 'CFD';
%input.F12 = 0.889304030368328;  % 0.889304030368328, 0.778608060736655
input = input_processing(input);
[t_FDM, t_FDM2, y_FDM, y_FDM2] = obtain_sol_rad(input);
[time_FDM, Temp_FDM, Tshelf_FDM, tm_FDM] = temperature_postprocessing(t_FDM, y_FDM, input);
[time_FDM2, itf_FDM, Tshelf_FDM2, tdry_FDM, Tb_FDM] = interface_postprocessingV2(t_FDM2, y_FDM2, input);
tdry_rad_sim1 = tdry_FDM;

end


%% Simplified Network
switch Sim_Multiple
case 'on'
    
input = input_def;
input.mode = 'CFD';
input.nx_vial = 10;
input.ny_vial = 10;
input.calmode = 'PreComputed';
input = input_processing(input);
input = get_ViewFactor_network(input);
F = input.F;
Nvial = input.Nvial;

% Simulation
tdry_rad_sim1 = zeros(Nvial ,1);
for i = 1:Nvial 
    input.F12 = F(i,end);
    input.rad = 'sim';
    input = input_processing(input);
    [t_FDM, t_FDM2, y_FDM, y_FDM2] = obtain_sol_rad(input);
    [time_FDM, Temp_FDM, Tshelf_FDM, tm_FDM] = temperature_postprocessing(t_FDM, y_FDM, input);
    [time_FDM2, itf_FDM, Tshelf_FDM2, tdry_FDM, Tb_FDM] = interface_postprocessingV2(t_FDM2, y_FDM2, input);
    tdry_rad_sim1(i) = tdry_FDM;
end

end


%% Radiation Network
switch Network
case 'on'
    
% Radiation input setup
input = input_def;
input.mode = 'CFD';
input.nx_vial = 10;
input.ny_vial = 10;
input.c = 0;
input.N = 3;
input.layout = [];
input.calmode = 'PreComputed';
input = input_processing(input);
input = get_ViewFactor_network(input);

% Simulation
[tdry_vial_network, sol_t, sol_T] = obtain_sol_rad_network(input);

end


