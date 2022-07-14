%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% The functions of the discretized radial displacement equation.
% Updated with backward-diff method instead of central diff method.
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       u       = radial displacment array
%       C       = concentrations at a time step 
%       Omega   = volumetric expansion factor 
%
% Returns:
%       F       = array of radial displacement equation functions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = u_funcs(N, R_oi, u, C, Omega)
    
    % Preallocation for speed
    F = zeros(1, N + 1);

    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial values
    R = 1:h:R_oi;

    % BC @ R = R_in
    F(1) = (1 + (u(1)/R(1))) - (1 + Omega * C(1))^(2/3); 
    
    % Eqs 2, 3, ... N+1
     for i = 2:N + 1
         F(i) = 1 + ((u(i) - u(i - 1))/(h)) - ((1 + Omega * C(i))/(1 + u(i)/R(i)));
     end 

end 