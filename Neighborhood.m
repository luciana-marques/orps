%===============================================================
% Demand Response Problem
% Title: Neighborhood Strucutre Generation
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Author: Luciana Sant'Ana Marques Arnoux and Isabella
% Date: Jun 19th, 2017 at 15:19
%===============================================================

function [loadsNew, appIndex] = Neighborhood(loads, nLoads)
% Input: 
    % loads: (1 x nLoads) structure with a solution of the loads
    % nLoads: (1 x 1) 
% Action:
    % 1- To interruptible:
    % Generate another feasible solution by choosing a consumer n and an 
    % interruptible appliance a at random and changing the status of the 
    % on/off variables at two time slots 
    
    % 2- To uninterruptible:
    % Generate another feasible solution by choosing a consumer n and an 
    % uninterruptible appliance a at random and changing the status of the 
    % on/off variables at two time slots    
% Output:
    % loadsNew: (1 x countLoads) structure with a new solution
    % appIndex: (1 x 1) index of the modified load

%Begin
    % New loads structure gets old loads structure
    loadsNew = loads;
    
    % Choose a load
    a = randi(nLoads,1,1);
    appIndex = a;
    
    % Take the solution vector
    sol = loads(a).solution;
    
    % Check which type of appliance
    
    % If interruptible load
    if loads(a).isUn == 0
        
        % Execute Nfirst
        
        % Sort the solution vector
        [B,I] = sort(sol);

        % Calculate how much positions are 0
        C = logical(B);
        nZeros = size(C(not([C])),2);

        % Calculate how much positions are 1
        nOnes = loads(a).beta - loads(a).alpha + 1 - nZeros;

        % If appliance has more possible time slots than its duration
        if (nZeros ~= 0) && (nOnes ~= 0)

            % Generate a random position for zeros
            posZeros = randi(nZeros,1,1);

            % Generate a random position for ones
            posOnes  = randi(nOnes,1,1);

            % Change variable status for posZeros and posOne at solution vector
            sol(I(posZeros)) = 1;
            sol(I(nZeros+posOnes)) = 0;

        end
            
    % If uninterruptible    
    else
        % Execute Nsecond
        
        % Take the duration vetor
        duration = loads(a).duration;

        % Option chose: Change only ONE position to deslocate the solution
        % Decision of direction
        if sol(1) == 1
            % Direction: right
            sol(duration+1) = 1;
            sol(1) = 0;
        else
            if sol(size(sol,2)) == 1
                % Direction: left
                sol(size(sol,2)-duration) = 1;
                sol(size(sol,2)) = 0;
            else
                aux = rand();
                if aux <= 0.5
                    % Direction right
                    i = 1;
                    while sol(i) ~= 1
                        i = i + 1;
                    end
                    sol(i + duration)= 1;
                    sol(i) = 0;    
                else
                    % Direction left
                    i = size(sol,2);
                    while sol(i) ~= 1
                        i = i - 1;
                    end
                    sol(i - duration)= 1;
                    sol(i) = 0; 
                end
            end
        end  
    end
    
    % New loads solution for appliance a gets new solution sol
    loadsNew(a).solution = sol;
    
end
    