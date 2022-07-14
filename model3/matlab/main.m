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
nu = 0.3;

% Settings 
N = 32; 
tspan = [0, 1]; 

% N = 32 on [0, 1]
t = readmatrix('data/times.csv');
C = readmatrix('data/C.csv');

% Plots results of u, P11, P22
plot_u(N, R_oi, tspan, t, C, Omega, nu)
plot_P11(N, R_oi, tspan, t, C, Omega, nu)
plot_P22(N, R_oi, tspan, t, C, Omega, nu)