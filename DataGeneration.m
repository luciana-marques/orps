%===============================================================
% Demand Response Problem
% Title: Data Generation
% Based on: Vasirani and Ossowski (2012)
% Institution: Federal University of Minas Gerais (UFMG)
% Department: Graduate Program in Electrical Engineering
% Author: Luciana Sant'Ana Marques Arnoux and Isabella
% Date: Jun 14th, 2017 at 11:54
%===============================================================

function [countLoads, loads, w, pi, pc, b, S] = DataGeneration(N, last, delta, seed)
% Input: 
    % N: number of consumers
    % last: time horizon
    % delta: size of slot in hours
    % seed: if -1, no seed; 
% Action:
    % Generate a random instance based on Vasirani and Ossowski (2012)
% Output:
    % countLoads: (1 x N) vector of how many shiftable loads each consumer
                % has
    % loads: (1 x countLoads) structure with all data needed to solve
                % the scheduling problem
    % w: (N x last) vector with base load of each consumer
    % pi: (1 x last) vector with prices (TOU tariffs)
    % pc: peak pricing
    % b: (1 x last) vector with IBR
    
    % Initialize seed (if any)
    if seed == -1
        rng('shuffle')
    else
        rng('default');
        rng(seed);
    end

    % Sets of loads
    % ---------------------------------
    % All of them are presented in the article. Not every consumer 
    % have all of them. This will be decided by probability

    % Types of shiftable interruptible loads
    IS = {'HEAT8' 'HEAT9' 'HEAT10' 'HEAT11' 'HEAT12' 'HEAT13' 'HEAT14'...
          'HEAT15' 'HEAT16' 'HEAT17' 'HEAT18' 'HEAT19' 'HEAT20' 'HEAT21'...
          'HEAT22' 'AC' 'EVN' 'EVM'};

    % Types of shiftable uninterruptible loads
    AS = {'WASHING' 'DRYER' 'DISHWASHER'};

    % Types of shiftable loads
    S = [AS,IS];

    % Types of basic (non-shiftable) loads
    NS = {'WATER' 'LIGHTING' 'KITCHEN' 'FRIDGE' 'FREEZER' 'OVEN'...
          'MICROWAVE' 'TV' 'DESKTOP' 'LAPTOT'};
    % ----

    % Price of energy/capacity in the day-ahead market
    % ---------------------------------
    % We consider an hypothetical case (need to improve it after). We define
    % a base price and 6 periods in a day, each of them having one base price
    % multiplier

    basePrice = .16;
    priceMultiplier = [1 1.1 1.3 1.2 2 1.2];
    hourInterval = [0 6 11 13 17 22 24];

    % Calculate price considering the time slot size
    pi = zeros(1,last);
    j = 1;
    
    for i = 1:last 
        
      % Calculate hour
      aux = i*delta;
      
      % Calculate price
      pi(1,i) = basePrice*priceMultiplier(1,j);
      
      % Verify if it needs to change period
      if (aux >= hourInterval(1,j+1)) 
          j = j+1;
      end
    end

    % Capacity price (slightly bigger than the article)
    pc = N*0.07; 
    
    % Inclining block rate price
    ibrMultiplier = 1.5;
    b = ibrMultiplier*pi;   
    % ----

    % Base Load
    % ---------------------------------
    % Calculate the load for each equipment and consumer and aggregate them
    % Based entirely on article

    % Total base load
    w = zeros(N,last);

    % Base load of each consumer and equipment
    auxW = zeros(N,size(NS,2),last);

    % Probability of having each equipment
    prob = [0.2 1 0.53 1 0.23 0.77 0.9 1 0.52 0.41];

    % Equipment power
    power = [1.2 .176 .0165 .08 .07 1.2 1.3 .01 .1 .02];

    % For each consumer
    for n = 1:N

      % Electric Water Heater
      i = 1;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 
        % Initial time slot
        j = randi([1 6]);
  
        while (j <= last)
            % Until 18 hrs, the water heater is turned on for 
            % 30 min every 2 hrs
            while (j/6 < 18)
                auxW(n,i,j) = power(1,i);
                auxW(n,i,j+1) = power(1,i);
                auxW(n,i,j+2) = power(1,i);
                j = j+12;
            end
            % After 18 hrs, the water heater is turned on for 
            % 30 min every 1 hr
            while (j/6 <= 24)
                auxW(n,i,j) = power(1,i);
                if ((j+1)/6 <= 24) auxW(n,i,j+1) = power(1,i); end
                if ((j+2)/6 <= 24)auxW(n,i,j+2) = power(1,i); end
                j = j+6;
            end
        end

      end % Otherwise consumer n doen't have water heater (zeros)

      % Lighting
      i = 2;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 

          % Average from spanish document http://www.ree.es/sites/default/files/downloadable/atlas_indel_ree.pdf
          % p. 64: lighting represents 32% of power load at 22hrs -> 0,176 kW 
          mu = power(1,i);
          % Variance equals 5% of average
          sigma = sqrt(mu*.05);
          % Consumption is normally ditributed in every time slot
          aux = max(normrnd(mu,sigma,1,last),0);
          auxW(n,i,:) = aux;

      end % Otherwise consumer n doen't have lighting (zeros)

      % Kitchen
      i = 3;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 
          
          % Average from spanish document 
          % http://www.ree.es/sites/default/files/downloadable/atlas_indel_ree.pdf
          % p. 64: lighting represents 32% of power load at 22hrs -> 0,176 kW 
          mu = power(1,i);
          % Variance equals 5% of average
          sigma = sqrt(mu*.05);
          % Consumption is normally ditributed in every time slot
          % Until 6 am consumer doesn't use kitchen
          aux = max(normrnd(mu,sigma,1,last-1/delta*6),0);
          auxW(n,i,1/delta*6+1:last) = aux;

      end % Otherwise consumer n doen't have kitchen (zeros)

      % Fridge
      i = 4;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 

          aux = repmat(power(1,i),1,last);
          auxW(n,i,:) = aux;

      end % Otherwise consumer n doen't have fridge (zeros)

      % Freezer
      i = 5;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 

          aux = repmat(power(1,i),1,last);
          auxW(n,i,:) = aux;

      end % Otherwise consumer n doen't have freezer (zeros)

      % Oven
      i = 6;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 

          % Oven is used 2 times a week; so probability of using it at 
          % a day is 2/7

          % In a day
          p1 = rand; 
          if (p1 <= 2/7)

            % At lunch time is 0.8
            p2 = rand;
            if (p2 <= .8)
                % Run between 20 min and 1hr - number of times slots
                slots = randi([20/(delta*60) 60/(delta*60)]);
                % Start time slot: 14 +- 1hr
                startSlot = randi([13*(1/delta) 15*(1/delta)]);
                % Power consumption
                for j = startSlot:startSlot+slots-1
                    auxW(n,i,j) = power(1,i);
                end
            end

            % At dinner time is 0.2
            p2 = rand;
            if (p2 <= .2)
                % Run between 20 min and 1hr - number of times slots
                slots = randi([20/(delta*60) 60/(delta*60)]);
                % Start time slot: 21 +- 1hr
                startSlot = randi([20*(1/delta) 22*(1/delta)]);
                % Power consumption
                for j = startSlot:startSlot+slots-1
                    auxW(n,i,j) = power(1,i);
                end
            end

          end % Otherwise not used at this day

      end % Otherwise consumer n doen't have oven (zeros)

      % Microwave
      i = 7;
      p = rand;

      % If consumer has the equipment
      if (p <= prob(i)) 

          % Possible time intervals they use it (for 10 minutes)
          startSlot = [randi([8*(1/delta) 10*(1/delta)]),
                        randi([10*(1/delta) 12*(1/delta)]),
                        randi([14*(1/delta) 16*(1/delta)]),
                        randi([19*(1/delta) 21*(1/delta)])];

          probMicrowave = [.12 .2 .25 .43];

          for j = 1:4  
              p1 = rand;
              if (p1 <= probMicrowave(j))
                  auxW(n,i,startSlot(j)) = power(1,i);
              end
          end

      end % Otherwise consumer n doen't have microwave (zeros)

      % TV, Desktop and computer
      for i = 8:10

          p = rand;

          % If consumer has the equipment
          if (p <= prob(i)) 

              % It is used twice a day: 14 +- 1hr and 20 +- 1hrs
              startSlot = [randi([13*(1/delta) 14*(1/delta)]),
                           randi([19*(1/delta) 21*(1/delta)])];

              % Duration between 1 and 3 hrs
              duration = [randi([1*(1/delta) 3*(1/delta)]),
                           randi([1*(1/delta) 3*(1/delta)])];

              % Power consumption
              for j = startSlot(1):startSlot(1)+duration(1)-1
                auxW(n,i,j) = power(1,i);
              end

              % Power consumption
              for j = startSlot(2):startSlot(2)+duration(2)-1
                auxW(n,i,j) = power(1,i);
              end


          end % Otherwise consumer n doen't have TV or desktop or laptop (zeros)

      end % TV, desktop and computer

    end % all consumers analyzed

    w(:,:) = sum(auxW,2);

    % Shiftable loads
    % ---------------------------------
    % Calculate the preferred times and durations for each equipment and 
    % consumer
    % Based entirely on article

    % Probability of having each equipment
    prob = [.93 .28 .53 .41 .49 .1];
    
    % Starting time of loads
    ts = zeros(N,size(S,2));

    % Ending time of loads
    tf = zeros(N,size(S,2));

    % Duration of loads
    d = zeros(N,size(S,2));

    % Power rate of loads
    L = zeros(N,size(S,2));

    % Equipment power
    power = [.19 1.24 .66 1 1.5 1.92];
    
    % Structure of loads
    loads = struct('n',{},'type',{},'alpha',{},'beta',{},'power',{},...
                    'duration',{},'isUn',{},'solution',{});

    % Count how many loads each consumer have
    countLoads = zeros(1,N);
    
    % Counter
    count = 0;
                
    % For each consumer
    for n = 1:N

        % Uninterruptible loads

        % Washing machine
        % Used 3 times a week, at morning (11+-1hr to 15+-1hr)
        % with prob 78% or at night (19+-1hr to 23+-1hr).
        % Durantion draw from ~U[1,2] hours

        i = 1;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))

            % Used 3 times a week
            p1 = rand;
            if (p1 <= 3/7)
                
                % Count loads to be shcedule and of consumer n
                count = count + 1;
                countLoads(1,n) = countLoads(1,n) + 1;

                % At morning or at nigth
                p2 = rand;

                if (p2 <= .78) % at morning

                    % Starting time
                    js = randi([10*1/delta 12*1/delta]);

                    % Ending time
                    jf = randi([14*1/delta 16*1/delta]);

                else % at night

                    % Starting time
                    js = randi([18*1/delta 20*1/delta]);

                    % Ending time
                    jf = randi([22*1/delta 24*1/delta]);

                end

                % Consumer
                loads(count).n = n;
                
                % Type
                loads(count).type = 'WASHING';
                
                % Starting time
                loads(count).alpha = js;
                ts(n,i) = js;

                % Ending time
                loads(count).beta = jf;
                tf(n,i) = jf;

                % Duration (1 to 2 hours)
                duration = randi([1*1/delta 2*1/delta]);
                loads(count).duration = duration;
                d(n,i) = duration;
                
                % Load
                loads(count).power = power(i);
                L(n,i) = power(i);
                
                % It is uninterruptible?
                loads(count).isUn = true;
                
                % Initialize solution vector 
                loads(count).solution = zeros(1,jf-js+1);

            end % Otherwise not used this day

        end % Otherwise consumer n doesn't have washing machine (zeros)


        % Dryer
        % Used 3 times a week, between 17+-1hr and 21+-1hr
        % Duration between 1 and 1.5 hours

        i = 2;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))

            % Used 3 times a week
            p1 = rand;
            
            if (p1 <= 3/7)
                
                % Count loads to be shcedule and of consumer n
                count = count + 1;
                countLoads(1,n) = countLoads(1,n) + 1;
                
                % Consumer
                loads(count).n = n;
                
                % Type
                loads(count).type = 'DRYER';
                
                % Starting time
                js = randi([16*1/delta 18*1/delta]);
                loads(count).alpha = js;
                ts(n,i) = js;

                % Ending time
                jf = randi([20*1/delta 22*1/delta]);
                loads(count).beta = jf;
                tf(n,i) = jf;

                % Duration (1 to 2 hours)
                duration = randi([1*1/delta 1.5*1/delta]);
                loads(count).duration = duration;
                d(n,i) = duration;
                
                % Load
                loads(count).power = power(i);
                L(n,i) = power(i);
                
                % It is uninterruptible?
                loads(count).isUn = true;
                
                % Initialize solution vector 
                loads(count).solution = zeros(1,jf-js+1);

            end % Otherwise not used this day

        end % Otherwise consumer n doesn't have dryer (zeros)

        % Dishwasher
        % Used 4 times a week, starting at 15+-1hr (50%) or at 19+-1hr
        % Finishing at 19+-1hr or 23+-1hr
        % Duration between 1 and 2 hours

        i = 3;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))

            % Used 4 times a week
            p1 = rand;
            
            if (p1 <= 4/7)
                
                % Count loads to be shcedule and of consumer n
                count = count + 1;
                countLoads(1,n) = countLoads(1,n) + 1;

                % At morning or at afternoon
                p2 = rand;

                if (p2 <= .50) % at afternoon

                    % Starting time
                    js = randi([14*1/delta 16*1/delta]);

                    % Ending time
                    jf = randi([18*1/delta 20*1/delta]);

                else % at night

                    % Starting time
                    js = randi([18*1/delta 20*1/delta]);

                    % Ending time
                    jf = randi([22*1/delta 24*1/delta]);

                end
                
                % Consumer
                loads(count).n = n;
                
                % Type
                loads(count).type = 'DISHWASHER';
                
                % Starting time
                loads(count).alpha = js;
                ts(n,i) = js;

                % Ending time
                loads(count).beta = jf;
                tf(n,i) = jf;

                % Duration (1 to 2 hours)
                duration = randi([1*1/delta 2*1/delta]);
                loads(count).duration = duration;
                d(n,i) = duration;
                
                % Load
                loads(count).power = power(i);
                L(n,i) = power(i);
                
                % It is uninterruptible?
                loads(count).isUn = true;
                
                % Initialize solution vector 
                loads(count).solution = zeros(1,jf-js+1);
  
            end % Otherwise not used this day

        end % Otherwise consumer n doesn't have dishwasher (zeros)

        % Interruptible loads

        % Heating
        % must be on for 10 min in every hour between 8am and 8pm
        % must be on for 30 min in every hour between 8pm and 11pm
        % we will consider the homes have only one heater (in the article
        % they can have between 1 and 3)

        % Position of load heat in vector of all loads
        i = 4;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))

            % For each hour
            for tau = 8:22
                
                % Count loads to be shcedule and of consumer n
                count = count + 1;
                countLoads(1,n) = countLoads(1,n) + 1;
                
                % Consumer
                loads(count).n = n;
                
                % Type
                loads(count).type = 'HEATING';

                % Starting time
                js = 1/delta*tau + 1;
                loads(count).alpha = js;
                ts(n,tau-i) = js;

                % Ending time
                jf = 1/delta*tau + 1/delta;
                loads(count).beta = jf;
                tf(n,tau-i) = jf;

                % Duration
                if (tau <= 20) 
                    loads(count).duration = 1;
                    d(n,tau-i) = 1;
                else
                    loads(count).duration = 3;
                    d(n,tau-i) = 3;
                end

                % Power
                loads(count).power = power(i);
                L(n,tau-i) = power(i);
                
                % It is uninterruptible?
                loads(count).isUn = false;
                
                % Initialize solution vector 
                loads(count).solution = zeros(1,jf-js+1);
                
            end % all heat loads calculated

        end % Otherwise consumer n doen't heating (zeros)

        % AC
        % must be on between 1pm and 6 pm for a duration 
        % depending on the amount of energy consumed

        i = 5;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))
            
            % Count loads to be shcedule and of consumer n
            count = count + 1;
            countLoads(1,n) = countLoads(1,n) + 1;
            
            % Consumer
            loads(count).n = n;
             
            % Type
            loads(count).type = 'AC';

            % Starting time
            js = 1/delta*13 + 1;
            loads(count).alpha = js;
            ts(n,i+14) = js;

            % Ending time
            jf = 1/delta*18 + 1/delta;
            loads(count).beta = jf;
            tf(n,i+14) = jf;

            % Duration
            % Amount of energy for AC (from 1.6 to 5.6 kWh uniformaly
            % distributed)
            energy = 1.6+(5.6-1.6).*rand;
            duration = round(energy/(power(i)*delta));
            loads(count).duration = duration;
            d(n,i+14) = duration;
            
            % Power
            loads(count).power = power(i);
            L(n,i+14) = power(i);
                
            % It is uninterruptible?
            loads(count).isUn = false;
                
            % Initialize solution vector 
            loads(count).solution = zeros(1,jf-js+1);

        end % Otherwise consumer n doesn't have AC (zeros)

        % EV
        % Owner arrives at 19+-1hr and has to get out to work at 8+-1hr
        % So 2 shiftable loads

        i = 6;
        p = rand;

        % If consumer has the equipment
        if (p <= prob(i))
            
            % Batterie size (Nissan Leaf) = 24kWh
            bat = 24;

            % SOC of batterie when owner arrives
            soc = .3+(.8-.3).*rand;

            % For night period (when owner arrives)
            
            % Count loads to be shcedule and of consumer n
            count = count + 1;
            countLoads(1,n) = countLoads(1,n) + 1;
            
            % Consumer
            loads(count).n = n;
             
            % Type
            loads(count).type = 'EV';
            
            % Starting time
            js = randi([18*1/delta 20*1/delta]);
            loads(count).alpha = js;
            ts(n,i+14) = js;

            % Ending time
            jf = last;
            loads(count).beta = jf;
            tf(n,i+14) = jf;
            
            % k1 for duration calculation
            k1 = jf - js + 1;

            % Power
            loads(count).power = power(i);
            L(n,i+14) = power(i);
                
            % It is uninterruptible?
            loads(count).isUn = false;
                
            % Initialize solution vector 
            loads(count).solution = zeros(1,jf-js+1);

            % For morning period (before owner has to go)
            
            % Count loads to be shcedule and of consumer n
            count = count + 1;
            countLoads(1,n) = countLoads(1,n) + 1;
            
            % Consumer
            loads(count).n = n;
             
            % Type
            loads(count).type = 'EV';
            
            % Starting time
            js = 1;
            loads(count).alpha = js;
            ts(n,i+15) = js;

            % Ending time
            jf = randi([7*1/delta 9*1/delta]);
            loads(count).beta = jf;
            tf(n,i+15) = jf;

            % Duration
            
            % k2 for duration calculation
            k2 = jf;

            energy = bat*(1-soc);

            e1 = energy*k1/(k1+k2); %night charge
            e2 = energy*k2/(k1+k2); %morning charge

            duration = round(e1/(power(i)*delta)); %night duration
            loads(count-1).duration = duration;
            d(n,i+14) = duration;
            duration = round(e2/(power(i)*delta)); %morning duration
            loads(count).duration = duration;
            d(n,i+15) = duration;
            
            % Power
            loads(count).power = power(i);
            L(n,i+15) = power(i);
                
            % It is uninterruptible?
            loads(count).isUn = false;
                
            % Initialize solution vector 
            loads(count).solution = zeros(1,jf-js+1);

        end % Otherwise consumer n doen't have EV (zeros)

    end % all consumers analyzed
                     
end