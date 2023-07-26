% ==============================================================================
% This is a function for generating the vials.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
function [xc,yc,c_x,c_y] = generate_vials_network(nx_vial,ny_vial,n_vial,bx,by,phi,d,c,layout)

nc = 1;
xc = zeros(n_vial,1);
yc = zeros(n_vial,1);
    
% For all vials
if strcmp('CP',layout)
    for i = 1:ny_vial
        for j = 1:nx_vial
            xc(nc) = 0.5*(bx(2)-(nx_vial-1)*d*cos(pi/6)-d) + d/2 + (j-1)*(d)*cos(pi/6);
            yc(nc) = 0.5*(by(2)-(ny_vial)*d-d*sin(pi/6)) + d/2 + (i-1)*(d) + (1-rem(j,2))*d*sin(pi/6);
            nc = nc+1;
        end
    end
else
    for i = 1:ny_vial
        for j = 1:nx_vial
            xc(nc) = 0.5*(bx(2)-nx_vial*d-(nx_vial-1)*c) + d/2 + (j-1)*(c+d);
            yc(nc) = 0.5*(by(2)-ny_vial*d-(ny_vial-1)*c) + d/2 + (i-1)*(c+d);
            nc = nc+1;
        end
    end
end

c_x = cell(n_vial,1);
c_y = c_x;
for i = 1:n_vial
    c_x{i}= xc(i)+(d/2)*cos(phi);
    c_y{i}= yc(i)+(d/2)*sin(phi);
end


return