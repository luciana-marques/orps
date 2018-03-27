%===============================================================
% Load Scheduling Problem with Quadratic Cost
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Author: Luciana Sant'Ana Marques
% Date: Feb 20th, 2018 at 16:31
%===============================================================

% -----------------------------------------------------
% Parameters to be set

format long
last = 144;
delta = 1/6;
R = ones(1,last);

% Groups of instances
N = [3 10 20 50 100 200 1000 2000 5000 10000]; % number of consumers at each instance
                                               % if want to generate a new
                                               % set. To change other
                                               % parameters go to
                                               % GenerateData.m function
   
% Pricing model
isPP = true;                   % do not change it (only PP was implemented)

% Generate or load instance set
generateInstance = false;      % true if want to generate instance
                               % false if want to load existent instance
fromSeed = false;              % true if want to choose seed to generate
                               % false if randomly generated instance
InstanceSeed = 1;              % choose seed for instance generation (if applicable)
instanceName = 'TesteInstances2.mat'; % choose instance name to be loaded (if applicable)

% Optimization
solveGurobi = true;           % if want to solve each instance using Gurobi via AMPL
                               % (need license for both) 
modelName = '201802ModelPeakPricing.mod'; % optimization model
saveResultsGurobi = true;      % save results gurobi

% Heuristics parameters
wLocal = true;                 % with local search after each SA best solution improvement
localAfterSA = true;           % with local search after each SA complete
maxInitial = 30;               % number of SA runs for each instance
heuristicsRandom = true;       % if want to run the heuristics randomly or using a seed
HeuristicsSeed = 6;            % choose seed for heuristic procedure (if applicable)
solveHeuristics = false;        % true if want to solve SA for the group of instances
loadResults = false;           % true if want to load previous results
resultsName = 'solution-woLocal-smallInstances2-30x.mat'; % choose results .mat name to load
resultsSaveName = 'solution-TesteInstances2'; % choose results name to save
saveResultsSA = true;          % save heuristics results 
alphaIC = .05;                 % type I error for IC calculation
CIlevel = alphaIC/2;           % Type I Error
timeLim = 900;                % Time limit for SA and gurobi

% Print convergence plot
printConvergence = false;      % true if want to print convergence plot at each SA run

% -----------------------------------------------------
% Generate instances w or w/o seed or 
% Load an instance set

if generateInstance
    if fromSeed
        instance = InstanceGeneration(N,last,delta,InstanceSeed);
    else
        instance = InstanceGeneration(N,last,delta,-1);
    end
else
    load(instanceName);
    % Get N from the instance loaded
    N = zeros(1,size(instance,2));
    for i = 1:size(N,2)
        N(1,i) = size(instance(i).count,2);
    end
end

% -----------------------------------------------------
% Heuristics parameters and resolution

% Define initial temperature for each instance
if solveHeuristics
    
    % Initialize Heuristics
    if heuristicsRandom
        rng('shuffle');
    else
        rng('default');
        rng(HeuristicsSeed);
    end

    t0 = zeros(1,size(N,2));
    p = .95;                    % initial acceptance rate

    repN = 1000;

    for n = 1:size(N,2)

        % Get instance values
        countLoads = instance(n).count;
        loadsOr = instance(n).loads;
        w = instance(n).w;
        pi = instance(n).pi;
        pc = instance(n).pc;
        b = instance(n).b;
        S = instance(n).S;


        loads = RankingHeuristic(N(1,n),last,delta,loadsOr,sum(w),pi,countLoads);

        [totalCost,loadCurve] = TotalCostF(last,delta,loads,w,pi,pc,b,...
                        N(1,n)*ones(1,last),isPP);

        costRef = totalCost;
        diffMean = zeros(1,repN);

        for i = 1:repN

            [laux, appIndex] = Neighborhood(loads,size(loads,2)); 

            [totalCost,auxloadCurve] = UpdateCost(last, delta, laux(appIndex), ...
                            loads(appIndex), pi, pc, b, N(1,n)*ones(1,last), ...
                            isPP, loadCurve);

            diff = max(totalCost - costRef,0);
            diffMean(1,i) = diff;
        end

        t0(1,n) = (-mean(diffMean))/log(p);
    end

    % Levels for parameters (result from parameters tuning approach)
    alpha = .9;
    Mk = 30*round(sqrt(N));  

    % Record solutions SA
    bestValue = Inf(size(N,2),maxInitial);
    time = zeros(size(N,2),maxInitial);
    par = zeros(size(N,2),maxInitial);
    bestSolutions = struct('loads', {}, 'loadCurve', {}, 'fo', {}, 'par', {}); 
    CI = zeros(2,size(N,2));
    CI_time = zeros(2,size(N,2));
    CI_par = zeros(2,size(N,2));

    for j=1:size(N,2)

        % Get instance values
        countLoads = instance(j).count;
        loadsOr = instance(j).loads;
        w = instance(j).w;
        pi = instance(j).pi;
        pc = instance(j).pc;
        b = instance(j).b;
        S = instance(j).S;

        % Number of replications of the procedure
        for i = 1:maxInitial
            
            tic

            loads = RankingHeuristic(N(1,j),last,delta,loadsOr,...
                sum(w),pi,countLoads);

            % Calculate total cost using function
            [totalCost,loadCurve] = TotalCostF(last,delta,loads,...
                w,pi,pc,b,N(1,j)*R,isPP);

            % Simulated Annealing
            [optLoads, optTotalCost, costs, loadCurve] = SimulatedAnnealing(t0(1,j),alpha,...
                    Mk(1,j),last,delta,loads,pi,pc,b,N(1,j)*R,isPP,totalCost,...
                    wLocal,localAfterSA,loadCurve, timeLim);

            % Calculate PAR
            par(j,i) = max(loadCurve)/mean(loadCurve);

            % Record solution
            bestValue(j,i) = optTotalCost;
            time(j,i) = toc;

            % Print
            X = ['N = ',num2str(N(1,j)),' initial = ',num2str(i),' PAR = ',...
                num2str(par(j,i)),' FO = ',num2str(bestValue(j,i)),...
                ' TIME = ',num2str(time(j,i))];
            disp(X)

            if printConvergence

                % Initialize bestCosts
                bestCosts = zeros(1,3);
                
                % For converge plot
                bestCosts = [bestCosts; costs'];

                % Plot convergence
                figure
                bestCosts(~any(bestCosts,2),:) = [];
                plot(bestCosts,'LineWidth',1.5)
                legend('current', 'previous', 'optimal')
                xlabel('iteration')
                ylabel('total cost ($)')
                print('resultConvergence' ,'-depsc')

            end

            % Best solution update
            if optTotalCost <= min(bestValue(j,:))
                bestSolutions(j).loads = optLoads;
                bestSolutions(j).loadCurve = loadCurve;
                bestSolutions(j).fo = optTotalCost;
                bestSolutions(j).par = par(j,i);
            end

        end
        
        % Calculate CI for results

        % CI of results
        m = bestValue(j,:)';
        SEM = std(m)/sqrt(length(m));                       % Standard Error
        TS = tinv([CIlevel  1-CIlevel],length(m)-1);        % T-Score
        CI(:,j) = repmat(mean(m),2,1) + TS'*SEM;            % Confidence Intervals

        % CI of time
        m = time(j,:)';
        SEM = std(m)/sqrt(length(m));                       % Standard Error
        TS = tinv([CIlevel  1-CIlevel],length(m)-1);        % T-Score
        CI_time(:,j) = repmat(mean(m),2,1) + TS'*SEM;       % Confidence Intervals
        
        % CI of par
        m = par(j,:)';
        SEM = std(m)/sqrt(length(m));                       % Standard Error
        TS = tinv([CIlevel  1-CIlevel],length(m)-1);        % T-Score
        CI_par(:,j) = repmat(mean(m),2,1) + TS'*SEM;        % Confidence Intervals

        % Save results
        if saveResultsSA
           auxN = N(j);
           auxName = [resultsSaveName '-N=' num2str(N(j))];
           auxBestValue = bestValue(j,:);
           auxTime = time(j,:);
           auxPar = par(j,:);
           auxBestSolutions = bestSolutions(j);
           auxInstance = instance(j);
           auxCI = CI(:,j);
           auxCItime = CI_time(:,j);
           auxCIpar = CI_par(:,j);
           auxMk = Mk(j);
           auxt0 = t0(j);
           save(auxName,'alpha','auxMk','auxt0','auxN','heuristicsRandom',...
                'HeuristicsSeed','generateInstance','InstanceSeed',...
                'fromSeed','auxBestValue','auxTime','auxPar',...
                'auxBestSolutions','auxInstance','maxInitial','alphaIC',...
                'auxCI','auxCItime','auxCIpar')
        end

    end
    
    % Save results
    if saveResultsSA
        save(resultsSaveName,'alpha','Mk','t0','N','heuristicsRandom',...
                'HeuristicsSeed','generateInstance','InstanceSeed',...
                'fromSeed','bestValue','time','par','bestSolutions',...
                'instance','maxInitial','alphaIC','CI','CI_time','CI_par')
    end
  
end

% Load results from previous SA runs (for the instance loaded/generated)
if loadResults
    load(resultsName);
end

% -----------------------------------------------------
% Solve using Gurobi - AMPL and Gurobi licenses needed

% From optimization (centralized)
centSolutions = struct('time', {}, 'fo', {}, 'peak', {}, 'consLoad', {},...
                'consCost', {}, 'bb', {}, 'gap', {});

% From optimization (decentralized)
loadCurveDec = zeros(size(N,2),last);
parDec = zeros(size(N,2),1);
bestValuesDec = zeros(size(N,2),1);
timeDec = zeros(size(N,2),1);

if solveGurobi
    
    % Configure matlab-AMPL API
    basef = fileparts(which('main'));
    addpath(fullfile(basef, '/amplapi/matlab'));
    setUp
  
    % Write data frames from instance (fixed for every instances)
    
    % Types of shiftable interruptible loads
    IS = {'AC'; 'EVN'; 'EVM'};

    % Types of shiftable uninterruptible loads
    AS = {'WASHING'; 'DRYER'; 'DISHWASHER'};
   
    % Types of shiftable loads
    S = [AS; IS];
      
    ampl = AMPL;                                 % Create AMPL object
    ampl.read(modelName);                        % Read optimization model
    ampl.setOption('solver', 'gurobi');          % Choose Gurobi as solver
    delete 'log.txt'                             % Delete log to write new
    contlim = 0;
    ampl.setOption('gurobi_options', 'timelim=900 logfile=log.txt mipgap=1e-4');
    
    ampl.getSet('IS').setValues(IS);             % Assign set IS
    ampl.getSet('AS').setValues(AS);             % Assign set AS
    ampl.getParameter('last').setValues(last);   % Assign last
    ampl.getParameter('delta').setValues(delta); % Assign delta
    
    for j=1:size(N,2)
        
        % Get data from instance
        pi = instance(j).pi;
        w = instance(j).w;
        n = N(1,j); 
        C = transpose(1:n);
        pc = instance(j).pc;
        
        % Construct matrices ts, tf, d and L
        auxLoads = instance(j).loads;
        nLoads = size(auxLoads,2);
        
        ts = zeros(n,size(S,1));
        tf = zeros(n,size(S,1));
        d = zeros(n,size(S,1));
        L = zeros(n,size(S,1));
        
        % For each load
        cont = 1;
        
        for k = 1:n
            
            if cont > nLoads
                    break
            end
        
            while auxLoads(cont).n == k

                if strcmp(auxLoads(cont).type, 'WASHING')
                    ts(k,1) = auxLoads(cont).alpha;
                    tf(k,1) = auxLoads(cont).beta;
                    d(k,1) = auxLoads(cont).duration;
                    L(k,1) = auxLoads(cont).power;
                    cont = cont + 1;
                elseif  strcmp(auxLoads(cont).type, 'DRYER')
                    ts(k,2) = auxLoads(cont).alpha;
                    tf(k,2) = auxLoads(cont).beta;
                    d(k,2) = auxLoads(cont).duration;
                    L(k,2) = auxLoads(cont).power;
                    cont = cont + 1;
                elseif  strcmp(auxLoads(cont).type, 'DISHWASHER')
                    ts(k,3) = auxLoads(cont).alpha;
                    tf(k,3) = auxLoads(cont).beta;
                    d(k,3) = auxLoads(cont).duration;
                    L(k,3) = auxLoads(cont).power;
                    cont = cont + 1;
                elseif  strcmp(auxLoads(cont).type, 'AC')
                    ts(k,4) = auxLoads(cont).alpha;
                    tf(k,4) = auxLoads(cont).beta;
                    d(k,4) = auxLoads(cont).duration;
                    L(k,4) = auxLoads(cont).power;
                    cont = cont + 1;
                elseif  strcmp(auxLoads(cont).type, 'EV')
                    for i = 5:6
                        ts(k,i) = auxLoads(cont).alpha;
                        tf(k,i) = auxLoads(cont).beta;
                        d(k,i) = auxLoads(cont).duration;
                        L(k,i) = auxLoads(cont).power;
                        cont = cont + 1;
                    end
                end
                
                if cont > nLoads
                    break
                end
            end
            
        end
        
        ampl.getParameter('cons').setValues(n);  % Assign number of consumers
        ampl.getParameter('pi').setValues(pi);   % Assign prices
        ampl.getParameter('w').setValues(w);     % Assign base load
        ampl.getParameter('ts').setValues(ts);   % Assign initial time
        ampl.getParameter('tf').setValues(tf);   % Assign final time
        ampl.getParameter('d').setValues(d);     % Assign duration
        ampl.getParameter('L').setValues(L);     % Assign power
        c = ampl.getSet('C');                    % Pointer for coalition
        c.setValues(C);                          % Assign consumers set
        ampl.getParameter('pc').setValues(pc);   % Assign peak prices
        
        tic
        ampl.solve()                             % Solve
        centSolutions(j).time = toc;             % Get resolution time
        
        z = ampl.getObjective('coalition_value');
        centSolutions(j).fo = z.getValues.getColumnAsDoubles('val');
                
        if strcmp({z.result},'limit')
            contlim = contlim + 1;
            fid = fopen('log.txt', 'r');
            s = textscan(fid, '%s', 'delimiter', '\n');
            aux_indx = find(strcmp(s{1}, 'Time limit reached'));
            indx = aux_indx(contlim);
            phrase = split(s{1}{indx+1},',');
            sphrase = split(phrase{2});
            centSolutions(j).bb = str2double([sphrase{4}]);
            sphrase = split(phrase{3});
            sphrase = split(sphrase{3},'%');
            centSolutions(j).gap = str2double([sphrase{1}])/100;
            fclose(fid);
        else
            centSolutions(j).bb = centSolutions(j).fo;
            centSolutions(j).gap = 0;
        end
        
        y = ampl.getVariable('P').getValues;     % get load curve
        total_load = y.getColumnAsDoubles('val');
        centSolutions(j).consLoad = reshape(total_load,[last,n]);
        
        peak = ampl.getVariable('xc').getValues; % get peak
        centSolutions(j).peak = peak.getColumnAsDoubles('val'); 
        
        total_energy = delta*centSolutions(j).consLoad; 
        costs = pi*total_energy + max(centSolutions(j).consLoad)*pc;
        centSolutions(j).consCost = costs*centSolutions(j).fo/sum(costs); % proportionally individual cost

        % Solve the model for each consumer independently
        tic
        for k = 1:n
          
          c.setValues(k);                              % Assign only consumer k to coalition
          ampl.solve();                                % Solve
          y = ampl.getVariable('P');                   % Creat object to access variables
          y = y.getValues;                             % Get value of variable
          values = y.getColumnAsDoubles('val');        % Transforme variable in MATLAB matrix
          load = transpose(sum(values,2));             % Transpose vector to 1 x 144
          loadCurveDec(j,:) = loadCurveDec(j,:)...     % Uptade total load curve
              + load;
          z = ampl.getObjective('coalition_value');    % Pointer to objective
          z = z.getValues.getColumnAsDoubles('val');   % Get objective value
          bestValuesDec(j,1) = bestValuesDec(j,1) + z; % Sum z to total decentralized cost

        end
        
        % Get processing time for n consumers
        timeDec(j,1) = toc;
         
        % Sum w
        parDec(j,1) = max(loadCurveDec(j,:))/mean(loadCurveDec(j,:));
        
    end
    
    ampl.close()
    
    % Save results
    if saveResultsGurobi
        prompt = {'Enter name you want for the workspace with Gurobi results:'};
        dlg_title = 'Input';
        answer = inputdlg(prompt,dlg_title);
        save(answer{1},'parDec','timeDec','bestValuesDec','centSolutions',...
            'loadCurveDec')
    end
   
end

% -----------------------------------------------------
% Calculate power curve without load scheduling
% Consumption starts at ts for all schedulable appliances

woLoadCurve = zeros(size(N,2),last);
woPAR = zeros(size(N,2),1);

for n = 1:size(N,2)
    
    % Number of schedulable loads
    nLoads = size(instance(n).loads,2);
    
    for i = 1:nLoads 
        
        % Get parameters for controllable load i
        taux = instance(n).loads(i).alpha;
        daux = instance(n).loads(i).duration;
        Laux = instance(n).loads(i).power;
        
        % Update load curve considering that appliance starts
        % at time taux (alpha)
        for j = taux:(taux+daux-1)
            woLoadCurve(n,j) = woLoadCurve(n,j) + Laux;
        end
        
    end
    
    % Sum base load
    woLoadCurve(n,:) = woLoadCurve(n,:) + sum(instance(n).w);
    
    % Calculate PAR
    woPAR(n,1) = max(woLoadCurve(n,:))/mean(woLoadCurve(n,:));
end

% -----------------------------------------------------
% Plots and tables

y = 1:last;

for n = 1:size(N,2)
    
    figure
    plot(y,woLoadCurve(n,:)',':',y,loadCurveDec(n,:),'--',...
        y,auxBestSolutions.loadCurve,...
        'LineWidth',1.1)

    legend({'W/O Load Scheduling' 'With Decentralized Load Scheduling'...
        'With Centralized Load Scheduling SA' 'With Centralized Load Scheduling Gurobi'},'Location','northwest')
    xlabel('time slot')
    ylabel('load (kW)')
    axis([0 144 10 65])

    print(['resultPlot' num2str(50)],'-depsc')
    
end