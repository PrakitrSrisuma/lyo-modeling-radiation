% ==============================================================================
% This is a function for obtaining the solution to the radiation network
% technique.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
function [tdry_vial_network, sol_t, sol_T] = obtain_sol_rad_network(input)


% Extract input data
N = input.N;  % number of nodes
Nvial = input.Nvial;  % number of vials in the network
endtime = input.endtime_d;  % final time
dt = input.dt;  % data collection frequency from the ODE solver
tol = input.tol;

% ODE solver setup
Tini_d = [input.T0_d*ones(Nvial*N,1); input.s0_d*ones(Nvial,1)];  % initial condition, shelf temperaute at the last position
tspan = (0:dt:endtime)';  % define the time span
S_check = zeros(Nvial,1);

% Extract the relevant matrices 
A_FDM = obtain_Jacobian_FDM(input);
M_FDM = eye(length(Tini_d));
options = odeset('Event', @(t_FDM,y_FDM) event_sublimation_completes_network(t_FDM,y_FDM,input), 'Mass', M_FDM, 'RelTol', tol, 'AbsTol', tol);
sol_t = [];
sol_T = [];
tm_network = zeros(Nvial,1);

tic
while min(S_check)<1-10*tol
    [t_FDM,y_FDM] = ode15s(@(t_FDM,y_FDM) obtain_ODEs_rad_network(t_FDM,y_FDM,A_FDM,tm_network,input), tspan, Tini_d, options); 
    sol_t = [sol_t;t_FDM];
    sol_T = [sol_T;y_FDM];

    for i = 1:Nvial
        k = find(sol_T(:,(i-1)*N+1)>=input.Tm_d-10*tol,1);

        if isempty(k)
        else
            tm_network(i,1) = sol_t(k);
        end
    end

    S_check = y_FDM(end,end-Nvial+1:end);
    tspan = (t_FDM(end):dt:endtime)';  % define the time span
    Tini_d = y_FDM(end,:);

    if sol_t(end)>=endtime
        break
    end

end
toc

sol_t = input.time_dim(sol_t)/3600;

tdry_vial_network = zeros(Nvial,1);
for i = 1:Nvial
    tdry_vial_network(i) = sol_t(find(sol_T(:,end-Nvial+i)>=1-10*tol,1));
end

sol_T = input.temp_dim(sol_T);

return