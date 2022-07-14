%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3
% Mathematical Modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find roots of eigenvalue eq. 
function roots = find_eigs(r_io, end_point)
    % Eigenvalue eq. 
    fun = @(mu) bessely(1, mu).*besselj(0,mu * r_io) - besselj(1, mu).*bessely(0, mu*r_io);
    format long;

    % Iterates through [0,end_point] range 
    % and finds all roots of the eq. 
    vals = []; 
    for i = 1:0.1:end_point
        val = fzero(fun, i);
        if isnan(val)
           continue;
        end
        vals = [vals, val]; 
    end
    
    % Return the unique elements that are to a unqiue TOL. 
    roots = uniquetol(vals, 10^(-7));
end 

