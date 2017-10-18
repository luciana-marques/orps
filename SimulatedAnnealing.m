% ====================================
% 16 June , 2017
% UFMG - PPGEE
% Optimal Residential Power Scheduling in Smart Grid
% Using Simulated Annealing to MINIMIZE the cost.
% Luciana e Isabella
% ====================================

function [optLoads, optTotalCost, costs, optLoadCurve] = SimulatedAnnealing(t0, alpha,...
                Mk, last, delta, loads, pi, pc, ro, R, isPP,...
                totalCost, wLocal,localAfterSA,loadCurve)
% Input:
    % loadCurve: matrix of appliances of n consumers
    % costShift: cost of shiftables 
    % t0
    % alpha
    % Mk
% Action:
    % Construct a solution using Simulated Annealing to
    % minimize cost.
    % Shake loads >> Calculate de current cost/ Compare/ Keep the min Cost or
    % by randon probability. Keep the best cost in "optimal_cost".
% Output:
    % loadCurve: Last loadCurve related to optimal_cost ?
    % totalCost: = optimal_cost of the solution
    
    % SA Parameters
    k = 1; % annealing parameter
    delta2 = 10^(-10);
    K = round((log(delta2)-log(t0))/log(alpha)); % number of temperatures
    stop = 0; % flag for stop criteria
    temperature = t0; 
    
    % Starting with the Previous Cost
    previous_cost = totalCost; 
    optimal_cost = totalCost; % Keep the best solution
    current_cost = totalCost;
    
    previous_load = loads;
    optimal_load = loads;
    previous_load_curve = loadCurve;
    
    % Counter of iterations without solution improvement
    iterations = 0;
    
    % Number of Temperatures without improvement
    numIt = 5;
    
    % Test Vector
    costs = zeros(3,K*Mk);
    
    % Load curve
    optLoadCurve = loadCurve;
    
    % Iterations counter
    i = 1;
    
    % Time limit
    timeLim = 900;
    tic

    % None of the stoping criteria and k <= K
    while (stop ~= 1) && (k <= K)

        m = 1;
        iterations = iterations + 1;

        while m <= Mk
            
            % Guarda resultados
            costs(1,i) = current_cost;
            costs(2,i) = previous_cost;
            costs(3,i) = optimal_cost;
            i = i + 1;
            
            % Shake the solution
            [loads, appIndex] = Neighborhood(previous_load,size(previous_load,2));
            
            % Calculate the current cost: 
            [current_cost,load_curve] = UpdateCost(last, delta, loads(appIndex), ...
                        previous_load(appIndex), pi, pc, ro, R, isPP, ...
                        previous_load_curve);
                    
            % Comparision to keep the best solution and matrix loadCurve:
            if current_cost <= optimal_cost  
                
                optimal_cost = current_cost;
                optimal_load = loads;
                optLoadCurve = load_curve;
                
                % If SA with local search
                if wLocal
               
                    % Perform local search to try to improve solution                        
                    [localSloads, localScost, localLoad_curve] = LocalSearch(last, delta, loads, ...
                            pi, pc, ro, R, isPP, current_cost, previous_load_curve);

                    % If Local Search returns a strictly better solution,
                    % get it as optimal solution
                    if localScost < current_cost
                        optimal_cost = localScost;
                        optimal_load = localSloads;
                        optLoadCurve = localLoad_curve;
                    end
                end
            end
                        
            % Calculate the difference 
            diff = current_cost - previous_cost;
            
            if diff <= 0
                
                previous_load = loads;
                previous_cost = current_cost;
                previous_load_curve = load_curve;
                
                % If solution was improved, iterations = 0.
                iterations = 0;
                
            else
                if rand(1) < exp(-diff/(temperature))
                    previous_load = loads;
                    previous_cost = current_cost;
                    previous_load_curve = load_curve;
                end
                
            end
            
            % Calculate time
            time = toc;
            
            % If solution is not improved until "numIt" iterations, stop algorithm 
            if (iterations > numIt) || (time > timeLim)
                stop = 1;
                break
            end

            m = m + 1;
        end
        
        % Uptade parameter k and next temperature
        k = k + 1;
        temperature = alpha*temperature;

    end
    
    if (stop ~= 1) && (localAfterSA == true)
        [optimal_load, optimal_cost, optLoadCurve] = LocalSearch(last, delta, optimal_load, ...
                            pi, pc, ro, R, isPP, optimal_cost, optLoadCurve);
    end

    % Final solution:
    optTotalCost = optimal_cost;
    optLoads = optimal_load;
    
end