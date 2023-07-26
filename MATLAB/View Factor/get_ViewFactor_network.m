% ==============================================================================
% This is a function for calculating the view factor.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = get_ViewFactor_network(input)

%% Input
d = input.d*1e2;  % diameter (cm)
c = input.c*1e2;  % distance between vials (cm)
nx_vial = input.nx_vial;  % number of vials per column
ny_vial = input.ny_vial;  % number of vials per row
Nvial = input.Nvial;
Nsurf = input.Nsurf;
phi = 0:pi/200:2*pi;
bx = [0;input.Lc*1e2];
by = [0;input.Lc*1e2];
layout = input.layout;
calmode = input.calmode;
Area = input.Area;

%% View Factor Calculation
switch calmode
% View factors for generic cases    
case 'General'

    % Create the vials
    [~,~,c_x,c_y] = generate_vials_network(nx_vial,ny_vial,Nvial,bx,by,phi,d,c,layout);
    
    % Intercept
    n_theta = 100;
    
    % View factor calculation
    n_count = 3;
    F = zeros(Nsurf,Nsurf,n_count);
    
    tic
    for i = 1:n_count
        [F(:,:,i),~,~] = cal_ViewFactor_network(c_x, c_y, bx, by, phi, Nvial, n_theta);
    end
    toc
    F = mean(F,3);
    
    % View factors for the wall
    for i = 1:Nsurf-1
        F(Nsurf,i) = F(i,Nsurf)*Area(i)/Area(Nsurf);
    end
    F(Nsurf,Nsurf) = 1-sum(F(Nsurf,1:end-1));

% View factors for one vial    
case '1Vial'
    F = zeros(Nsurf);
    F(1,1) = 0;  F(1,2) = 1;  
    F(2,1) = F(1,2)*input.A1/input.A2; F(2,2) = 1-F(2,1);

% View factors for two vials
case '2Vials'
    F = zeros(Nsurf);
    Y = 1 + input.c/input.d;
    F11 = (1/pi)*(sqrt(Y^2-1) + asin(1/Y) - Y);
    F12 = 1 - F11;

    F(1,1) = 0;  F(1,2) = F11;  F(1,3) = F12;  
    F(2,1) = F11;  F(2,2) = 0;  F(2,3) = F12;  
    F(3,1) = F(1,3)*input.A1/input.A2;  F(3,2) = F(2,3)*input.A1/input.A2;
    F(3,3) = 1 - F(3,1) - F(3,2);

case '3Vials'
    % View factors for three vials
    F = zeros(Nsurf);
    Y = 1 + input.c/input.d;
    F11 = (1/pi)*(sqrt(Y^2-1) + asin(1/Y) - Y);
    F12 = 1 - F11;
    F(1,1) = 0;  F(1,2) = F11;  F(1,3) = 0;  F(1,4) = F12;
    F(2,1) = F11;  F(2,2) = 0;  F(2,3) = F11;  F(2,4) = 1-2*F11;
    F(3,1) = 0;  F(3,2) = F11;  F(3,3) = 0;  F(3,4) = F12;
    F(4,1) = F(1,4)*input.A1/input.A2;  F(4,2) = F(2,4)*input.A1/input.A2;
    F(4,3) = F(3,4)*input.A1/input.A2;  F(4,4) = 1 - F(4,1) - F(4,2) - F(4,3);

case 'PreComputed'
    if strcmp('CP',layout)
        filename = ['ViewFactor_',num2str(nx_vial),'By',num2str(ny_vial),'_d',num2str(d),'_c',num2str(c),'CP'];
    else
        filename = ['ViewFactor_',num2str(nx_vial),'By',num2str(ny_vial),'_d',num2str(d),'_c',num2str(c)]; 
    end
    load(['Saved Data\ViewFactor\' filename])

case 'ChamberSize'
    filename = ['ViewFactor_',num2str(nx_vial),'By',num2str(ny_vial),'_d',num2str(d),'_c',num2str(c),'_side',num2str(rad.L)];
    load(['Saved Data\ChamberSize\' filename])
end

%% Export
input.F = F;
outputs = input;

return 

