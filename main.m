%===============================================================
% Load Scheduling Problem
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Course: Network Optimization
% Author: Luciana Sant'Ana Marques Arnoux and Isabella 
% Date: Jun 14th, 2017 at 11:54
%===============================================================

% -----------------------------------------------------
% Parameters to be set

format long
last = 144;
delta = 1/6;
R = ones(1,last);

% Groups of instances
N = [3]; % number of consumers at each instance
                                               % if want to generate a new
                                               % set. To change other
                                               % parameters go to
                                               % GenerateData.m function
   
% Pricing model
isPP = true;                   % do not change it (only PP was implemented)

% Generate or load instance set
generateInstance = true;       % true if want to generate instance
                               % false if want to load existent instance
fromSeed = false;              % true if want to choose seed to generate
                               % false if randomly generated instance
InstanceSeed = 1;              % choose seed for instance generation (if applicable)
instanceName = 'PSCC2018-Instances.mat'; % choose instance name to be loaded (if applicable)

% Optimization
solveGurobi = false;            % if want to solve each instance using Gurobi via AMPL
                               % (need license for both) 
modelName = '201705ModelPeakPricing.md'; % optimization model
saveResultsGurobi = true;      % save results gurobi

% Heuristics parameters
wLocal = true;                 % with local search after each SA best solution improvement
localAfterSA = true;           % with local search after each SA complete
maxInitial = 1;                % number of SA runs for each instance
heuristicsRandom = true;       % if want to run the heuristics randomly or using a seed
HeuristicsSeed = 3;            % choose seed for heuristic procedure (if applicable)
solveHeuristics = true;        % true if want to solve SA for the group of instances
loadResults = false;           % true if want to load previous results
resultsName = '';              % choose results .mat name to load
saveResultsSA = true;          % save heuristics results 
alphaIC = .05;                 % type I error for IC calculation

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

    % Initialize Heuristics
    if heuristicsRandom
        rng('shuffle');
    else
        HeuristicsSeed = 3; 
        rng('default');
        rng(HeuristicsSeed);
    end

    % Record solutions SA
    bestValue = Inf(size(N,2),maxInitial);
    time = zeros(size(N,2),maxInitial);
    par = zeros(size(N,2),maxInitial);
    bestSolutions = struct('loads', {}, 'loadCurve', {}, 'fo', {}, 'par', {}); 

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
                    Mk(1,j),last,delta,loads,pi,pc,b,N(1,j)*R,isPP,totalCost,wLocal,localAfterSA,loadCurve);

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

    end
   
    % Calculate CI for results

    % CI of results
    m = bestValue';
    CIlevel = alphaIC/2;                                % Type I Error  
    SEM = std(m)/sqrt(length(m));                       % Standard Error
    TS = tinv([CIlevel  1-CIlevel],length(m)-1);        % T-Score
    CI = repmat(mean(m),2,1) + TS'*SEM;                 % Confidence Intervals
     
    % CI of time
    m = time';
    SEM = std(m)/sqrt(length(m));                       % Standard Error
    TS = tinv([CIlevel  1-CIlevel],length(m)-1);        % T-Score
    CI_time = repmat(mean(m),2,1) + TS'*SEM;            % Confidence Intervals
    
    % Save results
    if saveResultsSA
        prompt = {'Enter name you want for the workspace with SA results:'};
        dlg_title = 'Input';
        answer = inputdlg(prompt,dlg_title);
        save(answer{1},'alpha','Mk','t0','N','heuristicsRandom',...
            'HeuristicsSeed','generateInstance','InstanceSeed',...
            'fromSeed','bestValue','time','par','bestSolutions',...
            'instance','maxInitial','alphaIC','CI','CI_time')
    end
end

% Load results from previous SA runs (for the instance loaded/generated)
if loadResults
    load(resultsName);
end

% -----------------------------------------------------
% Solve using Gurobi - AMPL and Gurobi licenses needed

% From optimization (centralized)
solutionGurobi = zeros(1,size(N,2));
timeGurobi = zeros(1,size(N,2));

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
    IS = {'HEAT8'; 'HEAT9'; 'HEAT10'; 'HEAT11'; 'HEAT12'; 'HEAT13'; 'HEAT14';...
          'HEAT15'; 'HEAT16'; 'HEAT17'; 'HEAT18'; 'HEAT19'; 'HEAT20'; 'HEAT21';...
          'HEAT22'; 'AC'; 'EVN'; 'EVM'};

    % Types of shiftable uninterruptible loads
    AS = {'WASHING'; 'DRYER'; 'DISHWASHER'};
   
    % Types of shiftable loads
    S = [AS; IS];

    % Types of basic (non-shiftable) loads
    NS = {'WATER'; 'LIGHTING'; 'KITCHEN'; 'FRIDGE'; 'FREEZER'; 'OVEN';...
          'MICROWAVE'; 'TV'; 'DESKTOP'; 'LAPTOT'};
      
    ampl = AMPL;                                 % Create AMPL object
    ampl.read(modelName);                        % Read optimization model
    ampl.setOption('solver', 'gurobi');          % Choose Gurobi as solver
    ampl.setOption('gurobi_options', 'timelim=900 logfile=logPP.txt');
    
    ampl.getSet('IS').setValues(IS);             % Assign set IS
    ampl.getSet('AS').setValues(AS);             % Assign set AS
    ampl.getSet('NS').setValues(NS);             % Assign set NS
    ampl.getParameter('last').setValues(last);   % Assign last
    ampl.getParameter('delta').setValues(delta); % Assign delta
    
    
    for j=1:size(N,2)
        
        % Get data from instance
        pi = instance(j).pi;
        w = instance(j).w;
        n = N(1,j); 
        C = transpose([1:n]);
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
                elseif  strcmp(auxLoads(cont).type, 'HEATING')
                    for i = 4:18
                        ts(k,i) = auxLoads(cont).alpha;
                        tf(k,i) = auxLoads(cont).beta;
                        d(k,i) = auxLoads(cont).duration;
                        L(k,i) = auxLoads(cont).power;
                        cont = cont + 1;
                    end
                elseif  strcmp(auxLoads(cont).type, 'AC')
                    ts(k,19) = auxLoads(cont).alpha;
                    tf(k,19) = auxLoads(cont).beta;
                    d(k,19) = auxLoads(cont).duration;
                    L(k,19) = auxLoads(cont).power;
                    cont = cont + 1;
                elseif  strcmp(auxLoads(cont).type, 'EV')
                    for i = 20:21
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
        timeGurobi(1,j) = toc;                   % Get resolution time
        
        z = ampl.getObjective('coalition_value').getValues;
        solutionGurobi(1,j) = z.getColumnAsDoubles('val');
        
        % Solve the model for each consumer independently
        tic
        for k = 1:n
        
            c.setValues(k);                              % Assign only consumer k to coalition
            ampl.solve();                                % Solve
            y = ampl.getVariable('y');                   % Creat object to access variables
            y = y.getValues;                             % Get value of variable
            values = y.getColumnAsDoubles('val');        % Transforme variable in MATLAB matrix
            values = reshape(values,...                  % Reshape matrix to 144
                [last,size(instance(1).S,2)]); 
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
        loadCurveDec(j,:) = loadCurveDec(j,:) + sum(instance(j).w);
        parDec(j,1) = max(loadCurveDec(j,:))/mean(loadCurveDec(j,:));
        
    end
    
    ampl.close()
    
    % Save results
    if saveResultsGurobi
        prompt = {'Enter name you want for the workspace with Gurobi results:'};
        dlg_title = 'Input';
        answer = inputdlg(prompt,dlg_title);
        save(answer{1},'parDec','timeDec','bestValuesDec','solutionGurobi',...
            'timeGurobi','loadCurveDec')
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
% Plots

y = 1:last;

for n = 1:size(N,2)
    
    figure
    plot(y,woLoadCurve(n,:)',':',y,loadCurveDec(n,:),'--',...
        y,bestSolutions(n).loadCurve,'LineWidth',1.1)

    legend({'W/O Load Scheduling' 'With Decentralized Load Scheduling'...
        'With Centralized Load Scheduling'},'Location','northwest')
    xlabel('time slot')
    ylabel('load (kW)')

    print(['resultPlot' num2str(N(1,n))],'-depsc')
    
end