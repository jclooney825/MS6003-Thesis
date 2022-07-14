%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3
% Mathematical Modeling MSc

% Returns the analytical solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical solution 
function [r,c] = analytic_sol(r_oi, num_terms, roots, t)
    % Radial points 
    r = linspace(1, r_oi, 1000); 
    
    % Initial condition 
    c = 1;

    % Bessel series summation 
    for i = 1:num_terms
        mu = roots(i);
        a_n = pi*(besselj(1, mu)*besselj(0, mu*r_oi))/( (besselj(1, mu))^2 - (besselj(0, mu*r_oi)^2));
        c = c - a_n*(besselj(0, mu*r).*bessely(1, mu) - besselj(1, mu).*bessely(0, mu*r))*exp(-(eigs(mu))^2*t);
    end 

end 