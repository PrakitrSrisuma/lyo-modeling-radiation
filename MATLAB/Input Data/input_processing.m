% ==============================================================================
% This is a function calculating all inputs.
% Lyophilization Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = input_processing(input)

% General Parameters/Variables
input.alpf = input.kf/(input.rhof*input.Cpf);  % thermal diffusivity (m2/s)
input.V = pi*input.d^2*input.L/4;  % volume of one vial (m3)
input.Tb0 = min(input.Tmax,input.Tb0);
input.tmax = max((input.Tmax-input.Tb0)/input.r, 0);  % time that the heater temperatura reaches its maximum value (s)

% Microwave
input.Qv = input.Q/input.V;  % power density
input.Hv1 = input.p1*input.Qv;
input.Hv2 = input.p2*input.Qv;
input.Hv3 = input.p3*input.Qv;  
input.Q_total = input.V*(input.Hv1 + input.Hv2 + input.Hv3);
assert(input.Q_total<input.Q,'Error, Microwave output exceeds the input power')

switch input.mode
case 'MFD'
input.hb = 0;
case 'CFD'
input.Hv1 = 0;
input.Hv2 = 0;
input.Hv3 = 0; 
end

% Sidewall radiation
if strcmp('sim',input.rad) || strcmp('est',input.rad)
    input.Nvial = 1;  % independent radiation
else
    input.Nvial = input.nx_vial*input.ny_vial;  % total number of vials
end
input.Nsurf = input.Nvial + 1;  % vials + 1 chamber wall
input.A1 = pi*input.d*input.L;  % surface area of one vial (m2)
input.A1_total = input.Nvial*input.A1;  % total surface area (m2)
input.V_total = input.Nvial*input.V;  % total volume (m3)
input.A2 = 6*input.Lc^2;  % surface area of the walls (m2)
if strcmp('est',input.rad)
    input.Rrad = input.Rrad;
else
    input.Rrad = (1-input.eps1)/(input.eps1*input.A1_total) + 1/(input.F12*input.A1_total) ...
     +(1-input.eps2)/(input.eps2*input.A2);
end

% Emissivity
input.eps = zeros(input.Nsurf,1);
for i = 1:input.Nsurf
    input.eps(i) = input.eps1;

    if i == input.Nsurf
        input.eps(i) = input.eps2;
    end
end

% Area
input.Area = zeros(input.Nsurf,1);
for i = 1:input.Nsurf
    input.Area(i) = input.A1;

    if i == input.Nsurf
        input.Area(i) = input.A2;
    end
end

% Functions
% non = nondimensionalize, dim = convert back to actual unit
input.temp_non = @(x) (x-input.Tb0)/(input.Tm-input.Tb0);
input.temp_dim = @(x) x*(input.Tm-input.Tb0)+input.Tb0;
input.temp_non2 = @(x) (x-input.Tmax)/(input.Tm-input.Tmax);
input.temp_dim2 = @(x) x*(input.Tm-input.Tmax)+input.Tmax;
input.pos_non = @(x) x/input.L;
input.pos_dim = @(x) x*input.L;
input.time_non = @(x) x*input.alpf/input.L^2;
input.time_dim = @(x) x*input.L^2/input.alpf;
input.irad_non = @(x) x*input.L^2/(input.kf*(input.Tm-input.Tb0));  % irradiation term
input.irad_dim = @(x) x*(input.kf*(input.Tm-input.Tb0))/input.L^2;
input.htc_non = @(x) x*input.L/input.kf;  % heat transfer coefficient or Biot number
input.htc_dim = @(x) x*input.kf/input.L;
input.shelf_non = @(x) x*input.L^2/(input.alpf*(input.Tm-input.Tb0));  % shelf heating rate
input.shelf_dim = @(x) x*(input.alpf*(input.Tm-input.Tb0))/(input.L^2);
input.sub_dim = @(x,y) input.L*(x+y*input.L)/(input.alpf*input.dHsub*(input.rhof-input.rhod));  % sublimation term
input.Tshelf = @(x) min([input.Tb0+input.r*x input.Tmax]);  %  shelf temperature (K)
input.Tshelf_d = @(x) min([input.shelf_non(input.r)*x input.temp_non(input.Tmax)]);  %  shelf temperature (dimensionless)
input.phi = @(x) input.L^2*x/(input.V_total*input.kf*(input.Tm-input.Tb0));
input.psi = @(x)input.L^2*x/(input.alpf*input.V_total*input.dHsub*(input.rhof-input.rhod));
input.phi_network = @(x) input.L^2*x/(input.V*input.kf*(input.Tm-input.Tb0));
input.psi_network = @(x)input.L^2*x/(input.alpf*input.V*input.dHsub*(input.rhof-input.rhod));

if strcmp('off',input.rad)
    input.phi = @(x) 0;
    input.psi = @(x) 0;
    input.phi_network = @(x) 0;
    input.psi_network = @(x) 0;
end

% Nondimensionalization
input.Tb0_d = input.temp_non(input.Tb0);
input.T0_d = input.temp_non(input.T0);
input.Tm_d = input.temp_non(input.Tm);
input.Ta_d = input.temp_non(input.Ta);
input.Tmax_d = input.temp_non(input.Tmax);
input.s0_d = input.pos_non(input.s0);
input.endtime_d = input.time_non(input.endtime);
input.tmax_d = input.time_non(input.tmax);
input.Bib = input.htc_non(input.hb);  % Biot number, nu in the manuscript
input.sigma = input.shelf_non(input.r);
input.lambda1 = input.irad_non(input.Hv1);
input.lambda2 = input.irad_non(input.Hv2);
input.lambda3 = input.irad_non(input.Hv3);

% Numerics and Discretization
input.dx = 1/(input.N-1);
input.x_d = (input.s0_d:input.dx:1)';
input.x = input.pos_dim(input.x_d);

% Export
outputs = input;

return