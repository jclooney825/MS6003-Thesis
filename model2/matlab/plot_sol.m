%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Plots the numerical solution on a given time range. 
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       tspan   = time range 
%       t       = times from solved equations
%       C       = solved concentrations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sol(N, R_oi, tspan, t, C)
    f1 =figure;
    
    % Step size 
    h = (R_oi - 1)/N; 
    
    % Radial values
    R = 1:h:R_oi;

    % Find indicies with t_val values 
    ind = []; 
    times = round(t, 2);
    t_vals = linspace((tspan(2)/4), tspan(2), 4); 
    
    for i = 1:1:4
        elements = find(times == t_vals(i));
        ind = [ind, elements(end)]; 
    end
    timestamps = t(ind); 

    % The concentration values at the selected times 
    y = C(ind,:); 
    
    % Plot solutions at various times 
    lines = {'--',':','-','-.'};
    hold on; 
    for i = 1:4
        %scatter(R, y(i,:),'black','filled', 'HandleVisibility', 'off');
        plot(R, y(i,:), 'DisplayName', strcat('t = ', num2str(times(i))),'LineWidth', 1.25, 'LineStyle', lines(i), 'Color', 'black');
    end
    xlim([1 R_oi])
    title('Numerical solutions of $C(r,t)$', 'Fontsize', 18, 'Interpreter','latex');
    xlabel('$r$', 'Interpreter','latex','Fontsize', 18);
    ylabel('$C$', 'Interpreter','latex', 'rotation', 0,'Fontsize', 18);
    str = strcat('$t = $', num2str(round(timestamps,2)));
    legend(str, 'Location','northwest', 'Fontsize', 12, 'Interpreter','latex');
end 