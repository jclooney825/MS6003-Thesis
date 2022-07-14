%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3
% Mathematical Modeling MSc

% Numerically solves the nondimensional circular diffusion equation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,c] = DiffusionSolver(r_oi, tspan, N)
    % Step size 
    h = (r_oi - 1) / N;

    % Radial iteration values 
    r_i = 1:h:r_oi;

    % A matrix construction 
    diags1 = [ones(N,1), -2*ones(N, 1), ones(N,1)];
    A = spdiags(diags1,-1:1, N, N);
    A(1,2) = 2;
    A = (1/h^2) * A; 

    % -1 and +1 off digonal values for matrix B
    off_diag_vals_min_1 = (1./circshift(r_i(1:end-1),-1))';
    off_diag_vals_plus_1 = (1./circshift(r_i(1:end-1),1))';

    % B matrix construction 
    diags2 = [-off_diag_vals_min_1 , zeros(N, 1), off_diag_vals_plus_1];
    B = spdiags(diags2,-1:1, N, N);
    B(1,2) = 0;
    B(end,end-1) = -(1/(r_i(end-1)));
    B = (1/(2*h)) * B;
   
    % epsilon vector 
    ep_vec = zeros(N, 1); 
    ep_vec(end) = (1/(h^2)) + (1/(2*h*r_i(end-1)));
        
    % Initial condition (t=0)
    c_init = zeros(N,1); 
    
    % Matrix equation derivatives 
    c_dot = @(t,c) (A + B) * c + ep_vec; 

    % ODE solver 
    options = odeset('RelTol',1e-13,'AbsTol',1e-20); 
    [t,c] = ode45(c_dot, tspan, c_init, options);
    c(:,end+1) = 1; 
end 


