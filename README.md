# ORPS
An optimization approach for residential power scheduling in Smart Grid.

Optimal Residential Power Scheduling is an approach to solve the optimization problem of a group of residential consumers, represented by an aggregator, in the day-ahead electricity market. Consumers have controllable and uncontrollable appliances to be scheduled centrally by the aggregator. They also have time preferences related to their appliances (constraints). We consider the peak pricing model for the decision making (objective function to be minimized).
The proposed approach consists of a problem specific construction algorithm, a Simulated Annealing (SA) refining method and a post-optimization local search. You will need Matlab to use/test it. You have to set some parameters at the beginning of the main.m code:

Groups of instances

•	N = [n_1 n_2 …. n_L]: if you want to generate L instances, you must define the number of consumers for each of them. For example, the vector N = [2 5 10] defines a group with three instances, with 2, 5 and 10 consumers each. If you want to load an existing instance, you can let it empty (N = [ ]). Data will be generated based on Spanish consumers from [1]. You can change consumer’s parameters in the “GenerateData.m” file;

Instances generation

•	generateInstance = true/false: set it true if you want to generate a new instance and false if you want to load an existing one;
•	fromSeed = true/false: set it true if you want to choose a seed to the instance generation and false if you want it randomly generated;
•	InstanceSeed = v: choose seed value v for instance generation (if “fromSeed = true”);
•	instanceName = 'name.mat': choose instance name to be loaded (if “generateInstance = false”);

Optimization parameters (need Gurobi and AMPL licenses, and AMPL API for Matlab)

•	solveGurobi = true/false: set it true if you want to solve each instance using Gurobi via AMPL and false otherwise; 
•	modelName = '201802ModelPeakPricing.md': this is the integer programming model proposed for exact optimization. You can change it if you want, but you will have to do modifications in the code (because it constructs the parametrization inside Matlab, via AMPL API);
•	saveResultsGurobi = true/false: set it true if you want to save the results of the exact optimization model. A dialog box will appear after the optimization to ask the name of the file you want to save;

Heuristics parameters

•	wLocal = true/false: set it true if you want to do the Local Search each time the solution is improved in the Simulated Annealing procedure and false otherwise;
•	localAfterSA = true/false: set it true if you want to do the Local Search after the SA (only in the best solution);
•	maxInitial = m: set the number of SA runs m for each instance;
•	heuristicsRandom = true/false: set it true if you want to run the heuristics randomly and false if you want to use a seed;
•	HeuristicsSeed = h: choose a seed h for heuristic procedure (if “heuristicsRandom = false”);
•	solveHeuristics = true/false: set it true if you want to solve the problem with the proposed heuristic approach and false otherwise;
•	loadResults = true/false: set it true if you want to load previous results and false otherwise;
•	resultsName = results.mat': choose “results .mat” name to load (if “loadResults = true”);
•	saveResultsSA = true/false: set it true if you want to save the SA results. A dialog box will appear after the SA optimization to ask the name of the file you want to save;
•	alphaIC = ic: select the type I error ic for confident interval calculation;

Others

•	printConvergence = true/false: set it true if you want to print convergence plot at each SA run and false otherwise.

[1] M.  Vasirani  and  S.  Ossowski,  “A  collaborative  model  for  participatory load management in the smart grid,” in Workshop on AI Problems and Approaches for Intelligent Environments, 2012, p. 21.
   
 
