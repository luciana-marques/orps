#=====================================================================
# Demand Response Problem
# Title: Load Scheduling Model v2 (w/o y[n,a,t])
# File: .mod
# Version: Binary Variables without considering consumers that do not
# participate in coalition
# Pricing model: peak pricing
# Institution: Federal University of Minas Gerais (UFMG)
# Department: Graduate Program in Electrical Engineering
# Prof: Wadaed Uturbey
# Author: Luciana Sant'Ana Marques Arnoux
# Date: Fev 16th, 2018 at 13:32
#=====================================================================

###  SETS AND PARAMETERS  ###

reset;

param cons > 0 integer; # Number of consumers of aggregator

set N := 1..cons;	# Set of consumers

set C within N;		# Set of the coalition

set IS;			# Set of interruptible shiftable appliances

set AS;			# Set of atomic shiftable appliances 

set S := AS union IS;	# Set of all shiftable appliances  

param last > 0 integer; # Number of time intervals in a day

set T := 1..last;   	# Set of time intervals in a day

param delta >= 0; 	# Duration of a time slot (in hours)

###  MARKET PARAMETERS  ###

param pi {T} >= 0;	
			# Energy price at time t ($/kWh)

param pc >= 0;	
			# Capacity price ($/kW)

###  CONSUMERS PARAMETERS  ###

param w {N,T} >= 0;
			# Base load of all consumer at each 
			# time interval (in kW)

param ts {N, S} >= 0;
			# Minimum starting time for shiftable 
			# loads of consumer n (in hours)

param tf {N, S} >= 0;
			# Maximum ending time for shiftable 
			# loads of consumer n (in hours)

param d {N, S} >= 0;
			# Duration of the loads of  
			# consumer n (in hours)

param L {N,S} >= 0;
			# Power rate of the loads of  
			# consumers (in kW)

###  VARIABLES  ###

var x {C, S, T}, binary;
			# Scheduled load of shiftable appliance 
			# of consumer n in C at each time interval
			# (in kW). 1 if load is on at time t
			# 0 otherwise

var xc >= 0; 		# Peak power consumption of the group 
			# (in kW)

var P {C, T} >= 0;
			# Scheduled power curve (in kW). 

###  OBJECTIVE  ###	

minimize coalition_value: sum{t in T} (pi[t]*delta*(sum{n in C} (w[n,t] + 
								  sum{a in S} x[n,a,t]*L[n,a]))) + 
								  pc*xc; 
			# The value of the coalition C is the    
			# total cost of the contracted amout
			# of energy w[n,t] + y[n,a,t] + z[n,t]  
			# plus the cost of the capacity

###  CONSTRAINTS  ###

subject to check_prefer {n in C, a in S}: 
			tf[n,a] - ts[n,a] + 1 >= d[n,a];
			# The interval must be greater than the
			# time duration for each shiftable appliance
			# a of each consumer n

subject to cons_prefer {n in C, a in S}: 
			sum{t in T: t >= ts[n,a] and t <= tf[n,a]} x[n,a,t] 
			= d[n,a];
			# The appliance energy consumption may
			# be distributed within the preferred  
			# time intervals of consumers.
			# With binary variables, the sum of 
			# on time slots must be equal the 
			# duration.

subject to cons_out_prefer {n in C, a in S, t in T: t < ts[n,a] 
			or t > tf[n,a]}: x[n,a,t] = 0;
			# Out of these time intervals the consumption
			# of the appliance is zero, in other words,
			# the appliance is off  

subject to cons_atomic {n in C, a in AS, k in 2..last-1, 
			t in T: t >= ts[n,a] and t <= tf[n,a]-k}: 
			x[n,a,t] - x[n,a,t+1] + x[n,a,t+k] <= 1;
			# Constraint to atomic load to be
			# consecutive. If at two slots (t and
			# t+k) x are equal to 1, then there is
			# no slot in between that x is equal to 0.

subject to cons_capacity {t in T}: 
			xc >= sum{n in C} (w[n,t] + sum{a in S} x[n,a,t]*L[n,a]);
			# Maximum capacity needed to respect
			# consumers schedule

subject to load_conversion {t in T, n in C}: 
			P[n, t] = w[n,t] + sum{a in S} x[n,a,t]*L[n,a]; 
			# Easier to read after