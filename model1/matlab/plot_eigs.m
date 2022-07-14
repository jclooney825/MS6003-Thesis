%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3
% Mathematical Modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the roots to ensure they make sense. 
function plot_eigs(r_oi, end_point)
    % Find the roots
    roots = find_eigs(r_oi, end_point); 

    % Define eigenvalue equations
    mu = 0:0.001:end_point; 
    fun = @(mu) (besselj(0, mu*r_oi).*bessely(1, mu) - besselj(1, mu).*bessely(0, mu*r_oi));
    format long;
    
    % Plot functions and roots
    f3 = figure;
    hold on; 
    plot(mu, fun(mu), 'LineWidth', 1.25, 'Color', 'Black');
    plot(mu, 0*mu, 'LineWidth', 1.5, 'Color', 'blue');
    scatter(roots,0*roots, 'Color','red','MarkerFaceColor', 'r')

   % Plot settings
   xlim([1 mu(end)])
   grid on;
   title('Roots of Eigenvalue Equation', 'Fontsize', 18, 'Interpreter','latex');
   xlabel('$\mu$', 'Interpreter','latex','Fontsize', 18);
end 