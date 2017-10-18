%===============================================================
% Ranking Heuristic to Load Scheduling Problem
% Construction of an initial solution
% Based on: Wu et al. (2012) - with many modifications
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Jun 12th, 2017 at 18:43
%===============================================================

function [loads] = RankingHeuristic(N, last, delta, loads, w, pi, aN)
% Input: 
    % N: number of consumers
    % last: time horizon
    % delta: size of slot in hours
    % loads: structure with loads information
    % w: (1 x last) vector with base load of all consumers
    % pi: (1 x last) vector with prices (R$/kWh)
    % aN: (1 x N) vector with how many loads each consumer has
% Action:
    % Construct a random solution to scheduling loads problem trying to
    % minimize peak-to-average ratio (PAR)
% Output:
    % loads: structure with loads information

    % Ranking process
    costShift = 0;
    loadCurve = zeros(N,last);
    auxPi = pi;
    maxPi = max(pi) + 1;

    % Randomize loads
    tLoads = sum(aN);
    a = 1:tLoads;
    lRand = a(randperm(length(a)));

    % For each load
    for i = 1:tLoads
        
        % Consumer of this load
        n = loads(lRand(i)).n;
        
        % Auxiliary variables
        auxAlpha = loads(lRand(i)).alpha;
        auxBeta = loads(lRand(i)).beta;
        auxDurat = loads(lRand(i)).duration;
        auxPower = loads(lRand(i)).power;

        % Get vector pi between starting and ending time
        pInterval = auxPi(1,auxAlpha:auxBeta);
        pIntervalaux = pi(1,auxAlpha:auxBeta);

        % Sort pi
        [B,I] = sort(pInterval);

        % If interruptible load
        if loads(lRand(i)).isUn == 0
            % Fill solution vector of load with 1 until durantion
            % respecting sorted solution
            for j = 1:auxDurat
                loads(lRand(i)).solution(I(j)) = 1;
            end

        % If uninterruptible    
        else
            % Fill solution vector of load with 1 until durantion
            % in sequencial manner, after the firt less expensive 
            % period that respects the last starting time 
            % (beta - duration + 1)
            k = 1;

            % While the cheaper periods are greater than 
            % the last possible starting time for 
            % uninterruptible load
            while I(k) > (auxBeta - auxAlpha - auxDurat + 2)
                k = k + 1;
            end

            for j = I(k):(I(k)+auxDurat-1)
                loads(lRand(i)).solution(j) = 1;
            end
        end

        % Calculate cost (se eu mudar o vetor pi ao longo do tempo, criar
        % vetor auxiliar que guarda o vetor original)
        costShift = costShift + pIntervalaux*...
                (loads(lRand(i)).solution*auxPower*delta)';

        % Add load to total load curve of consumer n 
        loadCurve(n,auxAlpha:auxBeta) = loadCurve(n,auxAlpha:auxBeta) +...
                loads(lRand(i)).solution*auxPower;

        % Calculate total load curve of consumers
        loadCurveT = sum(loadCurve) + w;

        % Minimization of PAR

        % Find peak from load curve
        [B I] = sort(loadCurveT,'descend');

        % Get maxLoad from total load curve
        maxLoad = loadCurveT(I(1));

        % Reload pi in auxPi
        auxPi = pi;

        % For the positions where load == maxLoad
        contPP = 1;
        while loadCurveT(I(contPP)) == maxLoad
            auxPi(1,I(contPP)) = maxPi;
            contPP = contPP + 1;
        end

    end
end
