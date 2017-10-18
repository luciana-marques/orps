%===============================================================
% Calculate Total Cost of a Solution
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Jun 19th, 2017 at 16:47
%===============================================================

function [totalCost,totalLoad] = TotalCostF(last, delta, loads, w, pi, pc, ro, R, isPP)
% Input: 
    % last: number of time slots
    % delta: size of time slot
    % loads: structure of solution with loads information
    % w: (1 x last) vector with base load of all consumers
    % pi: (1 x last) vector with prices (R$/kWh)
    % pc: peak charge
    % ro: (1 x last) price multiplier for IBR
    % R: (1 x last) load threshold
    % isPP: true if peak pricing and false if inclining block rates
% Action:
    % Calculates total cost of solution
% Output:
    % totalCost: total cost of solution
    
    % Quantity of loads in structure
    nApp = size(loads,2);
    
    % Total load
    totalLoad = zeros(1,last);
    
    % Calculate total load
    for i = 1:nApp
    
        % Auxiliary variables
        auxAlpha = loads(i).alpha;
        auxBeta = loads(i).beta;
        auxPower = loads(i).power;
        
        totalLoad(auxAlpha:auxBeta) = totalLoad(auxAlpha:auxBeta) + ...
            loads(i).solution*auxPower;
    end
    
    totalLoad = totalLoad + sum(w);
       
    % Calculate cost if peak pricing
    if isPP
        
        % Find peak
        peak = max(totalLoad);
        
        % Total cost
        totalCost = delta*totalLoad*pi' + peak*pc;
        
    % Calculate cost if IBR
    else
        
        totalCost = 0;
        
        % For each time slot verify if load exceeds threshold R
        for i = 1:last
            
            if totalLoad(i) <= R(i)
                totalCost = totalCost + totalLoad(i)*pi(i);
            else
                totalCost = totalCost + R(i)*pi(i) + ...
                    (totalLoad(i)-R(i))*pi(i)*ro(i);
            end
        end
        
    end
    
end