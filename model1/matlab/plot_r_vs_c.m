%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 3

% Plots the numerical and analytic solutions on the time interval [0,1].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [times] = plot_r_vs_c(r_io, N, t, c)
    f1 =figure;
    h = (r_io - 1) / N; 

    % Radial iteration values 
    r_num = 1:h:r_io; % i = 1, 2, ... N + 1

    % Find indicies with t_val values 
    ind = []; 
    times = round(t, 3); 
    t_vals = 0.25:0.25:1; 

    for i = 1:1:4
        elements = find(times == t_vals(i));
        ind = [ind, elements(end)]; 
    end
    timestamps = t(ind); 

    % The concentration values at the selected times 
    y = c(ind,:);
    
    % Roots in range [0,1000]
    end_point = 1000;
    roots = find_eigs(r_io, end_point);
    
    % Number of terms
    num_terms = 100; 

    lines = {'--',':','-','-.'};
    hold on; 
    for i = 1:1:4
         scatter(r_num', y, 'black', 'filled', 'HandleVisibility','off')
         [r,c] = analytic_sol(r_io, num_terms, roots, t_vals(i));
         plot(r,c, 'DisplayName', strcat('t = ', num2str(times(i))),'LineWidth', 1.25, 'LineStyle', lines(i), 'Color', 'black');
    end
    xlim([1 r_io])
    title('Analytical and Numerical solutions of $c(r,t)$', 'Fontsize', 18, 'Interpreter','latex');
    xlabel('$r$', 'Interpreter','latex','Fontsize', 18);
    ylabel('$c $', 'Interpreter','latex', 'rotation', 0,'Fontsize', 18);
    str = strcat('$t = $', num2str(round(timestamps,2)));
    legend(str, 'Location','southeast', 'Fontsize', 12, 'Interpreter','latex');
end 