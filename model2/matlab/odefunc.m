%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Produces differential equations needed for the concentration solver. 
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       u       = evaluated radial displacement 
%       C       = evaluated concentration 
%       Omega   = volumetric expansion factor 
%
% Returns:
%       dC       = differential eqution functions of dicretized 
%                  concentration equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dC = odefunc(N, R_oi, u, C, Omega)
    % Preallocation for speed 
    dC = zeros(1, N + 1); 

    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial values
    R = 1:h:R_oi;

    % BC at R = 1
    dC(1) = ((1 + u(1)/R(1))^2)/((1 + Omega*C(1))^3) * ((2 * C(2) - 2 * C(1))/(h^2));
    format long; 
    for i = 2:N
        u_R = 1 + u(i)/R(i); 
        y = 1 + Omega * C(i); 
       
        u_ddr = (u(i + 1) - 2 * u(i) + u(i - 1))/(h^2);

        C_dr = (C(i + 1) - C(i - 1))/(2 * h); 
        C_ddr = (C(i + 1) - 2 * C(i) + C(i - 1))/(h^2); 
        
        dC(i) = ((u_R)^2)/(y^3) * (C_ddr - (Omega/(y)) * (C_dr)^2 - (2*u_R/y) * u_ddr * C_dr + (1/R(i)) * C_dr);
    end 
    
   % BC at R = R_oi
   C_ddr_bc = 2 *(C(N) - C(N + 1))/h^2; 
   u_ddr_bc = (1/h^2) * (u(N) - u(N + 1) + h * (1 + Omega * C(N + 1))^(1/3) - h);
   dC(N + 1) = (1/(1 + Omega * C(N + 1))^(5/3)) * (C_ddr_bc + (2/h) - (Omega/(1 + Omega * C(N + 1))) - (2/(1 + Omega * C(N + 1))^(1/3)) * u_ddr_bc + 1/R(N + 1));

end 