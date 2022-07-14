%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% James Clooney 
% Week 6
% Mathematical Modeling MSc
%
% -----------------------------------------------------------------------------
% Plots the numerical solution of radial displacement on time interval.
% -----------------------------------------------------------------------------
%
% Arguments:
%       N       = number of spatial partitions
%       R_oi    = outer radius 
%       tspan   = time range 
%       t       = times from solved equations
%       C       = solved concentrations
%       Omega   = volumetric expansion factor 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_u(N, R_oi, tspan, t, C, Omega)
    f2 =figure;
    
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
    
    u0 = zeros(1, N+1); 
    u_options = optimset('Display','off');
    u = @(C) fsolve(@(u) u_funcs(N, R_oi, u, C, Omega), u0, u_options);
    
    % Plot solutions at different times
    lines = {'--',':','-','-.'};
    hold on; 
    for i = 1:4
        plot(R, u(y(i,:)), 'DisplayName', strcat('t = ', num2str(times(i))),'LineWidth', 1.25, 'LineStyle', lines(i), 'Color', 'black');
    end

    % Plot settings
    xlim([bounds(1) bounds(2)])
    title('Numerical solutions of $u(r,t)$', 'Fontsize', 18, 'Interpreter','latex');
    xlabel('$r$', 'Interpreter','latex','Fontsize', 18);
    ylabel('$u_r$', 'Interpreter','latex', 'rotation', 0,'Fontsize', 18);
    str = strcat('$t = $', num2str(round(timestamps,2)));
    legend(str, 'Location','northwest', 'Fontsize', 12, 'Interpreter','latex');
end 