%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% The functions of the discretized radial displacement equation.
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       u       = radial displacment array
%       C       = concentrations at a time step 
%       Omega   = xolumetric expansion factor 
%       nu      = Poisson's ratio 
%
% Returns:
%       F       = array of radial displacement equation functions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = u_functions(N, R_oi, u, C, Omega, nu)
    format long
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
    
    u_R_BC1 = u(1)/R(1); 
    x = 1 + Omega * C(1); 
    
    % First ghost point   
    u_G1 = u(2) - 2*h * sqrt(c * x^(2/3) - (1/b) * (1 + (1 + u_R_BC1)^2)) + 2*h; 
    
    % u derixatives at BC1
    du_BC1 = (u(2) - u_G1) / (2*h);
    ddu_BC1 = (u(2) - 2*u(1) + u_G1) / (h^2);
    
    % Elastic stain components 
    E11_BC1 = (1/2) * ((1 + du_BC1)^2 * x^(-2/3) - 1);
    E22_BC1 = (1/2) * ((1 + u_R_BC1)^2 * x^(-2/3) - 1);
    E33_BC1 = (1/2) * (x^(-2/3) - 1);
    
    % Free surface condition 
    P22_BC1 = x^(1/3) * a * (E11_BC1 + b * E22_BC1 + E33_BC1) * (1 + du_BC1);
    
    % Derivatives of elastic stain components
    dE11_BC1 = ((1 + du_BC1) * x^(-2/3)) * ddu_BC1; 
    dE22_BC1 = ((1 + u_R_BC1) * x^(-2/3)) * ((1/R(1)) * du_BC1 - (u_R_BC1)^2); 
    dE33_BC1 = 0; 
    
    d1 = b * E11_BC1 + E22_BC1 + E33_BC1; 
    d2 = b * dE11_BC1+ dE22_BC1 + dE33_BC1;
    
    dP11_BC1 = (a/3) * x^(-2/3) * ((3 * x * d2 * (1 + du_BC1)) + (3 * x * d1 * ddu_BC1));
    
    F(1) = dP11_BC1 - (P22_BC1 / R(1));

    % Eqs for i = 2,3,...N
    for i = 2:N
        
        % u derixatixes
        du = (u(i + 1) - u(i - 1)) / (2*h);
        ddu = (u(i + 1) - 2*u(i) + u(i - 1)) / (h^2);
        
        % C derivative
        dC = (C(i + 1) - C(i - 1)) / (2*h);
        
        u_R = u(i)/R(i);
        y = 1 + Omega * C(i); 

        % Elastic stain components 
        E11 = (1/2) * ((1 + du)^2 * y^(-2/3) - 1); 
        E22 = (1/2) * ((1 + u_R)^2 * y^(-2/3) - 1); 
        E33 = (1/2) * (y^(-2/3) - 1);
    
        % Piola-Kirchhoff tensor components 
        P11 = y^(1/3) * a * (b * E11 + E22 + E33) * (1 + du);
        P22 = y^(1/3) * a * (E11 + b * E22 + E33) * (1 + u_R);   
       
        % Derivatives of elastic stain components
        dE11 = ((1 + du) * y^(-2/3)) * ddu - (Omega/3) * (1 + du)^2 * y^(-5/3) * dC; 
        dE22 = ((1 + u_R) * y^(-2/3)) * ((1/R(i)) * du - (u_R)^2) - (Omega/3) * (1 + u_R)^2 * y^(-5/3) * dC; 
        dE33 = -(Omega / 3) * y^(-5/3) * dC; 
        
        e1 = b * E11 + E22 + E33; 
        e2 = b * dE11 + dE22 + dE33;
        
        % Derixatixes of P_11
        dP_11 = (a/3) * y^(-2/3) * ((Omega * dC * e1 * (1 + du))  + (3 * y * e2 * (1 + du)) + (3 * y * e1 * ddu)); 
    
        % Dixergence eq
        F(i) = dP_11 + ((P11 - P22) / R(i));
    
    end   
    
    u_R_BC2 = u(N + 1) / R(N + 1); 
    z = 1 + Omega * C(N + 1);
    
    % Second ghost point
    u_G2 = u(N) + 2*h * sqrt(c * z^(2/3) - (1/b) * (1 + (1 + u_R_BC2)^2)) - 2*h;
    
    % u derivatives at BC2
    du_BC2 = (u_G2 - u(N)) / (2*h);
    ddu_BC2 = (u_G2 - 2*u(N + 1) + u(N)) / (h^2);
    
    % Elastic stain components 
    E11_BC2 = (1/2) * ((1 + du_BC2)^2 * z^(-2/3) - 1);
    E22_BC2 = (1/2) * ((1 + u_R_BC2)^2 * z^(-2/3) - 1);
    E33_BC2 = (1/2) * (z^(-2/3) - 1);
    
    P22_BC2 = z^(1/3) * a * (E11_BC2 + b * E22_BC2 + E33_BC2) * (1 + du_BC2); 
    
    dE11_BC2 = ((1 + du_BC2) * z^(-2/3)) * ddu_BC2 - (Omega/3) * (1 + du_BC2)^2 * z^(-5/3); 
    dE22_BC2 = ((1 + u_R_BC2) * z^(-2/3)) * ((1/R(N + 1)) * du_BC2 - (u_R_BC2)^2) - (Omega/3) * (1 + u_R_BC2)^2 * z^(-5/3); 
    dE33_BC2 = -(Omega / 3) * z^(-5/3); 
    
    f1 = b * E11_BC2 + E22_BC2 + E33_BC2; 
    f2 = b * dE11_BC2 + dE22_BC2 + dE33_BC2;
    
    % Derixatixes of P_11
    dP11_BC2 = (a/3) * z^(-2/3) * ((Omega * f1 * (1 + du_BC2))  + (3 * z * f2 * (1 + du_BC2)) + (3 * z * f1 * ddu_BC2)); 

    F(N + 1) = dP11_BC2 - (P22_BC2 / R(N + 1));
    
end