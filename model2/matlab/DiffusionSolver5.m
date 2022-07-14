%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Concentration equation solver.
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       tspan   = time range 
%       u0      = initial value of radial displacement 
%       C0      = initial value of concentration 
%       Omega   = volumetric expansion factor 
%
% Returns:
%       t       = time values
%       C       = concentration values 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, C] = DiffusionSolver5(N, R_oi, tspan, u0, C0, Omega)
    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial values
    R = 1:h:R_oi;
    
    % Radial displacement function 
    u_options = optimset('Display','off');
    u = @(C) fsolve(@(u) u_funcs(N, R_oi, u, C, Omega), u0, u_options);
    
    % Time derivative 
    C_dot = @(t,C) odefunc(N, R_oi, u(C), C, Omega)';
    
    % Solve eqs and return 
    options = odeset('RelTol',1e-12,'AbsTol',1e-12); 
    [t,C] = ode45(C_dot, tspan, C0, options);

end 