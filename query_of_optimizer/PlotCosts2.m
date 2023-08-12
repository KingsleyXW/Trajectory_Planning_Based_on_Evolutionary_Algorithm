% Cite as:
% Mostapha Kalami Heris, Multi-Objective PSO in MATLAB (URL: https://yarpiz.com/59/ypea121-mopso), Yarpiz, 2015.


function PlotCosts(pop)

    pop_costs = [pop];
    plot(pop(:, 2), pop_costs(:, 1), 'bo');       
    xlabel('Left Area'); %1^{st} Objective 
    ylabel('Energy Consumption'); %2^{nd} Objective
    grid on;
    
end