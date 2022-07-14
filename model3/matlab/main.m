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

% Concentration profiles from model 2
t = readmatrix('data/t64.csv');
C = readmatrix('data/C64.csv');

% Settings 
N = width(C) - 1;
tspan = [t(1), t(end)];

% Plots results of u, P11, P22
plot_u(N, R_oi, tspan, t, C, Omega, nu)
plot_P11(N, R_oi, tspan, t, C, Omega, nu)
plot_P22(N, R_oi, tspan, t, C, Omega, nu)