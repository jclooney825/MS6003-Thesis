%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 7
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Plots the third diagonal component of the Piola-Kirchhoff stress tensor. 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_P33(N, R_oi, tspan, t, C, Omega, nu)
    f2 = figure; 

    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial values
    R = 1:h:R_oi;

    % Find indicies with t_val values 
    ind = []; 
    times = round(t, 4);
    t_vals = linspace((tspan(2)/4), tspan(2), 4);
    
    % Constants
    a = nu / ((1 + nu)*(1 - 2 * nu));
    b = (1 - nu) / nu; 
    
    for i = 1:1:4
        elements = find(abs(times - t_vals(i)) < 0.001);
        ind = [ind, elements(end)];
    end
    timestamps = t(ind);

    % The concentration values at the selected times 
    y = C(ind,:); 
    
    u0 = zeros(1, N+1); 
    u_options = optimset('Display','off');
    u = @(C) fsolve(@(u) u_functions(N, R_oi, u, C, Omega, nu), u0, u_options);
   
    lines = {'--',':','-','-.'};
    hold on; 
    for i = 1:4  
        P33 = calc_P33(N, R_oi, u(y(i,:)), y(i,:), Omega, nu);
        plot(R, P33, 'DisplayName', strcat('t = ', num2str(times(i))),'LineWidth', 1.25, 'LineStyle', lines(i), 'Color', 'black');
    end

    % Plot settings
    xlim([1 R_oi])
    title('Numerical solutions of $P_{33}$', 'Fontsize', 18, 'Interpreter','latex');
    xlabel('$R$', 'Interpreter','latex','Fontsize', 18);
    ylh = ylabel('$P_{33}$', 'Interpreter','latex', 'rotation', 0,'Fontsize', 18);
    ylh.Position(1) = ylh.Position(1) - 0.025 * ylh.Position(1); 
    str = strcat('$t = $', num2str(round(timestamps,2)));
    legend(str, 'Location','northeast', 'Fontsize', 12, 'Interpreter','latex');
    print('-depsc2','-painters','/Users/james/Desktop/MS6003-Thesis/model3/figures/eps_figs/P33.eps'); 
end 