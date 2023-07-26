% ==============================================================================
% This is a function for constructing a system of ODEs with the radiation
% network technique.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================

function outputs = obtain_ODEs_rad_network(t,y,A_FDM,tm,input)

% Extract some important input data
N = input.N;
Nvial = input.Nvial;
dx = input.dx;
Bib = input.Bib;
lambda1 = input.lambda1;
lambda3 = input.lambda3;
phi = input.phi_network;
psi = input.psi_network;
hb = input.hb;
Tm = input.Tm;
SB = input.SB;
Tw = input.Tw;
A = input.Area;
F = input.F;
eps = input.eps;
temp_dim = input.temp_dim;
temp_non = input.temp_non;
Tshelf = input.Tshelf(input.time_dim(t));
Tshelf_d = input.temp_non(Tshelf);


% Extract data for each vial
T = cell(Nvial,1);
S = zeros(Nvial,1);

for i = 1:Nvial
    T{i} = temp_dim(y((i-1)*N+1:i*N));
    S(i) = y(end-Nvial+i);

    if T{i}(1) >= Tm
        Tb = input.Tm_d + lambda3*(t-tm(i));
        T{i}(1:end) = input.temp_dim(Tb);
    end
end


% Network representation
Nsurf = Nvial + 1;
Qrad = zeros(N,Nsurf);

for i = 1:N
    Z = zeros(Nsurf );
    RHS_Z = zeros(Nsurf ,1);
    for j = 1:Nsurf
        if j == Nsurf
            RHS_Z(j) = -eps(j)*SB*Tw^4; 
        else
            RHS_Z(j) = -eps(j)*SB*T{j}(i)^4;
        end

        for k = 1:Nsurf
            if k == j
                Z(j,k) = -1+(1-eps(j))*F(j,k);
            else
                Z(j,k) = (1-eps(j))*F(j,k);
            end
        end

    end
    J = Z\RHS_Z;
    J = round(J,10);
%     CN = cond(Z);

    for j = 1:Nsurf-1
        Qrad(i,j) = (eps(j)*A(j)/(1-eps(j)))*(SB*T{j}(i)^4-J(j));  % portion of heat goes to the product
    end
    Qrad(i,end) = (eps(end)*A(end)/(1-eps(end)))*(SB*Tw^4-J(end));
end


% Create a RHS vector for the heating step
RHS = cell(Nvial,1);
for i = 1:Nvial
    RHS{i} = lambda1*ones(N,1);
    for j = 1:N
        RHS{i}(j) = RHS{i}(j) - phi(Qrad(j,i));
    end
    RHS{i}(end) = RHS{i}(end) + 2*dx*Bib*Tshelf_d/dx^2;
end


% Systems of equations
dydt = zeros(N*Nvial,1);
dsdt = zeros(Nvial,1);

for i = 1:Nvial
    if T{i}(1) >= Tm
        Hb = hb*(Tshelf-T{i}(1));
        Hv = input.Hv2;
        kappa = input.sub_dim(Hb, Hv);
        dsdt(i) = kappa + psi(-Qrad(1,i));
        dydt((i-1)*N+1:i*N) = zeros(N,1);
    else
        dsdt(i) = 0;
        dydt((i-1)*N+1:i*N) = A_FDM*temp_non(T{i}) + RHS{i};  
    end
end

outputs = [dydt;dsdt];

return