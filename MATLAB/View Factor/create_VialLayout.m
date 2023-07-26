% ==============================================================================
% This is a function for creating a figure of the vial layout.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function create_VialLayout(input)

%% Input
d = input.d*1e2;  % diameter (cm)
c = input.c*1e2;  % distance between vials (cm)
nx_vial = input.nx_vial;  % number of vials per column
ny_vial = input.ny_vial;  % number of vials per row
Nvial = input.Nvial;
phi = 0:pi/200:2*pi;
bx = [0;input.Lc*1e2];
by = [0;input.Lc*1e2];
layout = input.layout;

% Create the vials
[xc,yc,c_x,c_y] = generate_vials_network(nx_vial,ny_vial,Nvial,bx,by,phi,d,c,layout);
    
% % Check the center position
% plot(xc,yc,'.')
% xlim([0 bx(2)])
% ylim([0 by(2)])
% xline(bx(2)/2)
% yline(by(2)/2)
% close

% Plot the vial layout
for i = 1:Nvial
    plot(c_x{i},c_y{i},'-b','linewidth',1)
    hold on
end
xlim([bx(1) bx(2)])
ylim([by(1) by(2)])
set(gca,'linewidth',2,'XTick',[], 'YTick', [])
axis equal

return