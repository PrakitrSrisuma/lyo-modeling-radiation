% ==============================================================================
% This is a function creating and defining all inputs.
% Lyophilization Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = get_input_data

% Properties and Operating Conditions
input.mode = 'CFD';  % choose between 'HFD', 'CFD', 'MFD'
input.rad = 'on';  % 'on' for radiation network, 'sim' for simplfied network, 'est' for estimating Rrad, 'off' for no radiation
input.rhod = 63;  % density of dried layer (kg/m3)
input.rhof = 917;  % density of the frozen region (kg/m3)
input.Cpf = 1967.8;  % heat capacity of the frozen region (J/kgK)
input.kf = 2.30;  % thermal conductivity of the frozen region (kg/m3)
input.dHsub = 2.84e6;  % enthalpy of sublimation (J/kg)
input.hb = 65;  % heat transfer coefficient at the bottom surface (w/m2K), default 65
input.d = 0.01;  % vial diameter (m)
input.c = 0.005;  % distance between vials (m)
input.L = 0.042;  % length/height of the frozen material (m)
input.Tmax = 281.85;  % maximum shelf temperature (K), default 281.85
input.Tb0 = 236.85;  % initial shelf temperature (K), default 236.85
input.T0 = 236.85;  % initial temperature (K), default 236.85
input.Tm = 256.15;  % sublimation temperature (K)
input.Ta = 298.15;  % ambient temperature (K)
input.PT = @(x) exp(-6139.9/x + 28.8912);  % vapor pressure from temperature
input.TP = @(x) 6139.9*(1/(-log(x)+28.8912));  % sublimation temperature from pressure
input.s0 = 0;  % length/height of the frozen area (m)
input.r = 1/60;  % rate of shelf temperature increase (K/s), default 1/60

% Microwave
input.Q = 85;
input.p1 = 3.73e-04;
input.p2 = 8.62e-03;
input.p3 = 2.5e-05;

% Radiation
input.nx_vial = 1;  % number of vials in x-direction
input.ny_vial = 1;  % number of vials in y-direction
input.eps1 = 0.8;  % emissivity of glass
input.eps2 = 0.3;  % emissivity of stainless steel
input.SB = 5.67e-8;  % Stefan-Boltzmann constant
input.Lc = 0.3;  % chamber side (m)
input.Tw = 293.15;  % wall temperature (K)
input.F12 = 1;  % view factor; vial to wall
input.F = [];
input.calmode = 'General';
input.layout = [];

% Numerics and Discretization
input.endtime = 1e5;  % final time (s) 
input.tol = 1e-6;  % tolerance
input.N = 20;  % number of nodes
input.dt = 0.05;

% Export
outputs = input;

return