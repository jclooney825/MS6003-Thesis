%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Main function that calls all other functions. 
% -----------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
R_oi = 2; 
Omega = 10^(-3);

% Partition number 
N = 64;       

% Time range 
tspan = [0, 1];

% Initial conditions 
C0 = zeros(1, N+1); 
u0 = zeros(1, N+1); 

% Solution 
[t, C] = DiffusionSolver5(N, R_oi, tspan, u0, C0, Omega);

% Plot solutions 
plot_sol(N, R_oi, tspan, t, C)
plot_u(N, R_oi, tspan, t, C, Omega)