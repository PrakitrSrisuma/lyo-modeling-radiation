function [objective, termination, direction] = event_sublimation_completes_network(t,T,input)

N = input.N;
Nvial = input.Nvial;
obj = zeros(2*Nvial,1);
term = ones(2*Nvial,1);
dir = ones(2*Nvial,1);

for i = 1:Nvial
    obj(i) = T(N*(i-1)+1) - input.Tm_d + input.tol;
end
for i = 1:Nvial
    obj(end-Nvial+i) = T(end-Nvial+i) - 1 + input.tol;
end

objective = obj;  % stop when reaching sublimation temperature
termination = term;  % terminate ode solvers 
direction = dir;  % both directions

end