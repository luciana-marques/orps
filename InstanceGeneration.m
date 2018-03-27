%===============================================================
% Load Scheduling Problem - Instance Generation
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Author: Luciana Sant'Ana Marques
% Date: Feb 20th, 2018 at 16:33
%===============================================================

function instance = InstanceGeneration(N, last, delta, seed)
% Input: 
    % N: vector with number of consumers
    % last: time horizon
    % delta: size of slot in hours
    % seed: if -1, no seed; 
% Action:
    % Generate a set of random instances (.dat and struct)
% Output:
    % instance: struct with the data

    instance = struct('count',{},'loads',{},'w',{},'pi',{},'pc',{},...
                    'b',{},'alpha', {}, 'beta', {}, 'gamma', {}, ...
                    'S',{}, 'totalLoad', {});

    for i = 1:size(N,2)

        [countLoads, loadsOr, w, pi, pc, b, S, alpha, beta, gamma] = ...
            DataGeneration(N(1,i),last, delta, seed);

        instance(i).count = countLoads;
        instance(i).loads = loadsOr;
        instance(i).w = w;
        instance(i).pi = pi;
        instance(i).pc = pc;
        instance(i).b = b;
        instance(i).alpha = alpha;
        instance(i).beta = beta;
        instance(i).gamma = gamma;
        instance(i).S = S;
        instance(i).totalLoad = zeros(1,last);

    end

end