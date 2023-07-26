% ==============================================================================
% This is a function that calculated the view factor via Monte Carlo simulation.
% Freeze-drying Problem
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
function [F,n_ray,n_strike] = cal_ViewFactor_network( c_x, c_y, bx, by, phi, n_vial, n_theta)

n_high = 1e6;
F = zeros(n_vial+1, n_vial+1);

for k = 1:n_vial
    n_ray = 0;
    n_strike = zeros(length(phi),n_vial+1);
    parfor i = 1:length(phi)-1
        count = zeros(n_vial+1,1);
        x_in = n_high*ones(n_vial+1,1);
        y_in = n_high*ones(n_vial+1,1);
        theta = 0 + (pi)*rand(n_theta,1);
        %theta = 0 + linspace(0,pi,n_theta);

        for j = 1:n_theta
            n_ray = n_ray + 1;
            f = @(x,i) (tan(theta(j)-pi/2 + phi(i)))*(x-c_x{k}(i))+c_y{k}(i); %#ok<*PFBNS> 
            x_ray = [c_x{k}(i);c_x{k}(i)+bx(2)*cos(theta(j)-pi/2 + phi(i))];
            y_ray = f(x_ray,i);
%             plot(x_ray,y_ray)
%             hold on
            
            % Border
            [x0,y0,~,~] = intersections(x_ray,y_ray,[bx(1),bx(1),bx(2),bx(2),bx(1)],[by(1),by(2),by(2),by(1),by(1)],0);
            if isempty(x0)
                x_in(end) = n_high;
                y_in(end) = n_high;
            else
                x_in(end) = x0;
                y_in(end) = y0;           
            end
    
            % Vial
            for m = 1:n_vial
                if m ~= k
                    [x0,y0,~,~] = intersections(x_ray,y_ray,c_x{m},c_y{m},0);
                    if isempty(x0)
                        x_in(m) = n_high;
                        y_in(m) = n_high;
                    elseif length(x0) == 1
                        x_in(m) = x0;
                        y_in(m) = y0; 
                    else
                        d1 = norm([x0(1)-c_x{k}(i),y0(1)-c_y{k}(i)],2);
                        d2 = norm([x0(2)-c_x{k}(i),y0(2)-c_y{k}(i)],2);
                        
                        if d1 < d2
                            x_in(m) = x0(1);
                            y_in(m) = y0(1);
                        else
                            x_in(m) = x0(2);
                            y_in(m) = y0(2);
                        end
                    end
                end
            end

            dis = sqrt((x_in-c_x{k}(i)).^2+(y_in-c_y{k}(i)).^2);
            [~,index] = min(dis);

            count(index) = count(index) + 1;

        end
        n_strike(i,:) = count;
    end
    tmp = sum(n_strike);
    F(k,:) = tmp/n_ray;
end


% F = n_strike/n_ray;

return