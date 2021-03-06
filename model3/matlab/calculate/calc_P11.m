%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Calculates P11. Uses a central difference approx. for u derivative. 
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       u       = radial displacment array
%       C       = concentrations at a time step 
%       Omega   = volumetric expansion factor 
%       nu      = Poisson's ratio 
%
% Returns:
%       P11     = calculates P11 of Piola-Kirchhoff stress tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P11 = calc_P11(N, R_oi, u, C, Omega, nu)
   
    % Preallocation for speed
    F = zeros(1, N + 1); 

    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial xalues
    R = 1:h:R_oi;
    
    % Constants
    a = nu / ((1 + nu) * (1 - 2*nu));
    b = (1 - nu) / nu;
    c = (1 + nu) / (1 - nu);
    
    % Free surface condition 
    P11(1) = 0;

    % Eqs for i = 2,3,...N
    for i = 2:N
        
        % u derixatixes
        du = (u(i + 1) - u(i - 1)) / (2*h);
       
        u_R = u(i)/R(i);
        y = 1 + Omega * C(i);

        % Elastic stain components 
        E11 = (1/2) * ((1 + du)^2 * y^(-2/3) - 1);
        E22 = (1/2) * ((1 + u_R)^2 * y^(-2/3) - 1);
        E33 = (1/2) * (y^(-2/3) - 1);
    
        % Piola-Kirchhoff tensor components 
        P11(i) = y^(1/3) * a * (b * E11 + E22 + E33) * (1 + du);
    
    end 
    % Free surface condition 
    P11(N + 1) = 0;

end 