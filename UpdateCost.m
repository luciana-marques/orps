%===============================================================
% Calculate Total Cost of a Solution Considering only the
% Modified load
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Set 18th, 2017 at 17:06
%===============================================================

function [totalCost,totalLoad] = UpdateCost(last, delta, loads, ...
                        previous_Load, pi, pc, ro, R, isPP, ...
                        previousTL)
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
    % totalLoad
    
    % Previous solution vector
    pSol = previous_Load.solution;
    
    % New solution vector
    nSol = loads.solution;
    
    % Previous total load
    totalLoad = previousTL;

    % Auxiliary variables
    auxAlpha = loads.alpha;
    auxBeta = loads.beta;
    auxPower = loads.power;
    
    % Subtract previous load of appliance changed
    totalLoad(auxAlpha:auxBeta) = totalLoad(auxAlpha:auxBeta) - ...
            pSol*auxPower;
        
    % Add new solution
     totalLoad(auxAlpha:auxBeta) = totalLoad(auxAlpha:auxBeta) + ...
            nSol*auxPower;
    
    % Calculate cost if peak pricing
    if isPP
        
        % Find peak
        nPeak = max(totalLoad);
        
        % Total cost (- previous peak + new peak)
        totalCost = delta*totalLoad*pi' + nPeak*pc;
        
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