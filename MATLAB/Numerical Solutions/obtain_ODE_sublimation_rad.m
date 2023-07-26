% ==============================================================================
% This is a function creating the ODE for the sublimation stage.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = obtain_ODE_sublimation_rad(t,y,tm,input)

% Parameters
hb = input.hb;
Tm = input.Tm;
Tshelf = input.Tshelf(input.time_dim(t));
Tshelf_d = input.temp_non(Tshelf);
psi = input.psi;
Rrad = input.Rrad;
SB = input.SB;
Tw = input.Tw;
temp_dim = input.temp_dim;

% Create the ODE for sublimation
Hb = hb*(Tshelf-input.temp_dim(y(end-1)));
Hv = input.Hv2;
kappa = input.sub_dim(Hb, Hv);
Qrad = SB*(Tw^4-temp_dim(y(end-1))^4)/Rrad;
dsdt = kappa + psi(Qrad);

% Product temperature
Tb = input.Tm_d + input.lambda3*(t-tm);

outputs = [dsdt; y(end-1)-Tb; y(end)-Tshelf_d];

return