%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 8
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Main function that calls all other functions. 
% -----------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
R_oi = 2; 
Omega = 10^(-3);
nu = 0.3;

% Partition number 
N = 32;       

% Time range 
tspan = [0, 1];

% Initial conditions 
C0 = zeros(1, N+1); 
u0 = zeros(1, N+1); 

% Concentration profiles from model 2
[t, C] = DiffusionSolver5(N, R_oi, tspan, u0, C0, Omega);

% Plots results of u, P11, P22, P33
plot_u(N, R_oi, tspan, t, C, Omega, nu)
plot_P11(N, R_oi, tspan, t, C, Omega, nu)
plot_P22(N, R_oi, tspan, t, C, Omega, nu)
plot_P33(N, R_oi, tspan, t, C, Omega, nu)