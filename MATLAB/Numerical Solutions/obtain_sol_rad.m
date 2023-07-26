% ==============================================================================
% This is a function creating the FDM solution.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function [output1, output2, output3, output4] = obtain_sol_rad(input)

% Parameters
N = input.N;  % number of nodes
endtime = input.endtime_d;  % final time
dt = 0.05;  % data collection frequency from the ODE solver
Tini_d = [input.T0_d*ones(N,1); input.Tb0_d];  % initial condition, shelf temperaute at the last position
tspan = (0:dt:endtime)';  % define the time span

% Extract the relevant matrices 
A_FDM = obtain_Jacobian_FDM(input);
M_FDM = obtain_MassMatrix_FDM(input);
M_sub = zeros(3);
M_sub(1,1) = 1;

% Solve the system of ODEs using FDM // stage 1: increasing temperature
options_FDM = odeset('Event', @(t_FDM,y_FDM) event_sublimation_starts(t_FDM,y_FDM,input),'Mass', M_FDM, 'RelTol', 1e-9, 'AbsTol', 1e-9);
[t_FDM,y_FDM] = ode15s(@(t_FDM,y_FDM) obtain_ODEs_heating_rad(t_FDM,y_FDM,A_FDM,input), tspan, Tini_d, options_FDM); 

% Solve the system of ODEs using FDM // stage 2: sublimation
s0_FDM = [input.s0_d; y_FDM(end,end); input.Tshelf_d(t_FDM(end))];
options_sub = odeset('Event', @(t_FDM2,y_FDM2) event_sublimation_completes(t_FDM2,y_FDM2,input),'Mass', M_sub, 'RelTol', 1e-9, 'AbsTol', 1e-9);
[t_FDM2,y_FDM2] = ode15s(@(t_FDM2,y_FDM2) obtain_ODE_sublimation_rad(t_FDM2,y_FDM2,t_FDM(end),input), (t_FDM(end):dt:2*endtime)', s0_FDM, options_sub);

% Export
output1 = t_FDM;
output2 = t_FDM2;
output3 = y_FDM;
output4 = y_FDM2;

return