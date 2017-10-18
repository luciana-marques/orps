%===============================================================
% Load Scheduling Problem - Instance Generation
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Jul 12th, 2017 at 11:53
%===============================================================

function instance = InstanceGeneration(N, last, delta, seed)
% Input: 
    % N: vector with number of consumers
    % last: time horizon
    % delta: size of slot in hours
    % seed: if -1, no seed; 
% Action:
    % Generate a set of random instances (struct)
% Output:
    % instance: struct with the data

    instance = struct('count',{},'loads',{},'w',{},'pi',{},'pc',{},...
                    'b',{},'S',{}, 'totalLoad', {});

    for i = 1:size(N,2)

        [countLoads, loadsOr, w, pi, pc, b, S] = DataGeneration(N(1,i),...
            last, delta, seed);

        instance(i).count = countLoads;
        instance(i).loads = loadsOr;
        instance(i).w = w;
        instance(i).pi = pi;
        instance(i).pc = pc;
        instance(i).b = b;
        instance(i).S = S;
        instance(i).totalLoad = zeros(1,last);

    end

end