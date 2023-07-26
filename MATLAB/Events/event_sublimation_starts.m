% ==============================================================================
% This is a function defining the event where sublimation starts.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function [objective, termination, direction] = event_sublimation_starts(t,T,input)

objective = T(1)-input.Tm_d;  % stop when reaching sublimation temperature
termination = 1;  % terminate ode solvers 
direction = 0;  % both directions

end