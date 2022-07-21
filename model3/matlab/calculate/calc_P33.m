%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 7
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Calculates P33. Uses a central difference approx. for u derivative. 
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
%       P33     = calculates P33 of Piola-Kirchhoff stress tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P33 = calc_P33(N, R_oi, u, C, Omega, nu)
   
    % Preallocation for speed
    P33 = zeros(1, N + 1); 

    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial xalues
    R = 1:h:R_oi;
    
    % Constants
    a = nu / ((1 + nu) * (1 - 2*nu));
    b = (1 - nu) / nu;
    c = (1 + nu) / (1 - nu);
      
    u_R_BC1 = u(1)/R(1); 
    x = 1 + Omega * C(1); 
    
    % First ghost point   
    u_G1 = u(2) - 2*h * sqrt(c * x^(2/3) - (1/b) * (1 + (1 + u_R_BC1)^2)) + 2*h; 
    
    % u derixatives at BC1
    du_BC1 = (u(2) - u_G1) / (2*h);
    
    % Elastic stain components 
    E11_BC1 = (1/2) * ((1 + du_BC1)^2 * x^(-2/3) - 1);
    E22_BC1 = (1/2) * ((1 + u_R_BC1)^2 * x^(-2/3) - 1);
    E33_BC1 = (1/2) * (x^(-2/3) - 1);
    
    % Free surface condition 
    P33(1) = x^(1/3) * a * (E11_BC1 + E22_BC1 + b * E33_BC1);

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
        P33(i) = y^(1/3) * a * (E11 + E22 + b * E33);
    
    end 
     
    u_R_BC2 = u(N + 1) / R(N + 1); 
    z = 1 + Omega * C(N + 1);
    
    % Second ghost point
    u_G2 = u(N) + 2*h * sqrt(c * z^(2/3) - (1/b) * (1 + (1 + u_R_BC2)^2)) - 2*h;
    
    % u derivatives at BC2
    du_BC2 = (u_G2 - u(N)) / (2*h);
    
    % Elastic stain components 
    E11_BC2 = (1/2) * ((1 + du_BC2)^2 * z^(-2/3) - 1);
    E22_BC2 = (1/2) * ((1 + u_R_BC2)^2 * z^(-2/3) - 1);
    E33_BC2 = (1/2) * (z^(-2/3) - 1);
    
    P33(N + 1) = z^(1/3) * a * (E11_BC2 + E22_BC2 + b * E33_BC2); 
    
end 