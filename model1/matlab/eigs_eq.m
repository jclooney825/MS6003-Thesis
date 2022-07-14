%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 2
% Mathematical Modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigs_eq()
    % Define eigenvalue equations
    mu = 0:0.01:10; 
    fun = @(mu, r_oi) bessely(1,mu).*besselj(0,mu * r_oi) - besselj(1,mu).*bessely(0,mu * r_oi) ;
    format long;
    
    % Plot function
    f3 = figure;
    hold on;
   
    % Plots the eig eq. for various r_oi values
    for i = 2:1:5
        plot(mu, fun(mu, i), 'LineWidth', 1.25, 'DisplayName', strcat('$r_{oi} =  $',num2str(i)))     
    end
    
    % Plot settings
    xlim([0.5 mu(end)])
    grid on;
    title('Eigenvalue Equation Plot for various $r_{oi}$ values', 'Fontsize', 18, 'Interpreter','latex');
    legend('Location','southeast', 'Fontsize', 15, 'Interpreter','latex');
    xlabel('$\mu$', 'Interpreter','latex','Fontsize', 18);
    print('-depsc2','-painters','/Users/james/Desktop/Final Project/Week3/figures/eigs.eps'); 
end 