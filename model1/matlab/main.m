%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3
% Mathematical Modeling MSc

% Solutions to the eigenvalue problem and the numerical and analytical
% solutions of the radial diffusion equation. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10;           % Number of partitions 
tspan = [0, 1];   % Time range 
r_oi = 2;         % Ratio of radii 

% Numerical solution 
[t,c] = DiffusionSolver(r_oi, tspan, N);

% Plot solution 
plot_r_vs_c(r_oi, N, t, c);

% Plot roots of eig eq.
plot_eigs(r_oi, 50);

% Plot eig eq for various r_io
eigs_eq(); 