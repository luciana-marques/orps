%===============================================================
% Load Scheduling Problem - Local Search
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Sep 8th, 2017 at 11:40
%===============================================================

function [loads, fbest, loadCurve] = LocalSearch(last, delta, loads, ...
                        pi, pc, ro, R, isPP, totalCost, loadCurve)
% Input:
    % 
% Action:
    % 
% Output:
    % 

    % Best cost 
    fbest = totalCost;
    nLoads = size(loads,2);
    
    for a = 1:nLoads
        
        sol = loads(a).solution;
        auxSol = sol;
        
        % Check which type of appliance
        % If interruptible load
        if loads(a).isUn == 0
            
            for i = 1:size(sol,2)
                for j = (i+1):size(sol,2)
                    sol = auxSol;
                    if sol(i) ~= sol(j)
                        % Take the solution vector
                        aux = sol(i);
                        sol(i) = sol(j);
                        sol(j) = aux;
                        auxLoads = loads;
                        auxLoads(a).solution = sol;
                        
                        [auxCost, auxTotalLoad] = UpdateCost(last, delta, auxLoads(a), ...
                                        loads(a), pi, pc, ro, R, isPP, ...
                                        loadCurve);
                    
                        if auxCost <= fbest
                            fbest = auxCost;
                            loads(a).solution = sol;
                            loadCurve = auxTotalLoad;
                        end
                    end
                end
            end
            
        % If uninterruptible load
        else
            
            d = loads(a).duration;
            for i = 1:d
                if sol(i) == 1
                    break
                end
            end
            
            posIn = i;
            posFin = i+d-1;
            
            % Go left
            for i = 1:(posIn-1)
                % Take the solution vector
                sol(posIn-i) = 1;
                sol(posFin-i+1) = 0;
                auxLoads = loads;
                auxLoads(a).solution = sol;
                
                [auxCost, auxTotalLoad] = UpdateCost(last, delta, auxLoads(a), ...
                                        loads(a), pi, pc, ro, R, isPP, ...
                                        loadCurve);
                            
                if auxCost <= fbest
                    fbest = auxCost;
                    loads(a).solution = sol;
                    loadCurve = auxTotalLoad;
                end
            end
            
            % Size of soluction vector
            auxSize = size(sol,2);
            
            % Go right
            sol = auxSol;
            for i = (posFin+1):auxSize
                sol(i) = 1;
                sol(i-d) = 0;
                auxLoads = loads;
                auxLoads(a).solution = sol;
                
                [auxCost, auxTotalLoad] = UpdateCost(last, delta, auxLoads(a), ...
                                        loads(a), pi, pc, ro, R, isPP, ...
                                        loadCurve);                      
                        
                if auxCost <= fbest
                    fbest = auxCost;
                    loads(a).solution = sol;
                    loadCurve = auxTotalLoad;
                end
            end
            
        end
    end
end