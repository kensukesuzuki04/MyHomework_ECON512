% Empirical method HW4
% Kensuke Suzuki
% Penn State
% October 20

clear all
delete HW4log.txt
diary('HW4log.txt')
diary on

disp('ECON512 HOMEWORK4: Ken Suzuki')

%% Questtion 1: Quasi-Monte Carlo method with Dart Throwing

% Number of random draw
numsim = 10000;

seed = 1534561;
rng(seed);

% Define function: pi_ind
% returns 1 if x^2 + y^2 <= 1 and 0 otherwise

% Generate random sequence for x and y using rand
%seq = rand(numsim,2);
seq = haltonseq(numsim,2);
x = seq(:,1);
y = seq(:,2);

% We now compute the sequence of values of indicator function using the
% random sequence generated above
pi_QMC_Q1 = pi_ind(x,y);

% Compute numerical integation 
display('Problem 1: Quasi-Monte Carlo method')
pi_Q1 = ((1-0)*(1-0))/numsim * 4 * sum(pi_QMC_Q1)

clear x y
%% Questtion 2: Newton-Cotes approach with Dart Throwing
% Here I use the midpoint rule to compute the integration

% define the width of partition h
h = (1-0)/numsim;

% Define vector of x and y (later filled)
x = zeros(numsim,1);
y = zeros(numsim,1);

% x_j = 0 + (j-1/2)h for j=1,...,numsim
for ind = 1:numsim
    x(ind,1) = 0 + (ind- 1/2)*h;
    y(ind,1) = 0 + (ind- 1/2)*h;
end

% Compute the sequence of values of indicator function for given x_j
% and sum over with weight h, which yields the approximation of integration
% over y (for iven x_j)
pi_NC_Q2 = ones(numsim,1);
for ind = 1:numsim
    x_1 = x(ind,1)*ones(numsim,1);
    pi_NC_x = pi_ind(x_1,y);
    pi_NC_Q2(ind,1) =  h * sum(pi_NC_x);
end

% Next we integrate over x by summing over with weight h
display('Problem 2: Newton-Cotes approach')
pi_Q2 = 4 * h * sum(pi_NC_Q2)


%% Questtion 3: Newton-Cotes approach: Pythagorean
% I use Halton sequence to generate random draws

% Define function: pi_root
% returns (1-x^2)^(1/2)

seed = 1534561;
rng(seed);
% Generate random sequence for x 
x = haltonseq(numsim,1);

% We now compute the sequence of values of indicator function using the
% random sequence generated above
pi_QMC_Q3 = pi_root(x);

% Compute numerical integation 
display('Problem 3: Pythagorean fomula with Quasi-Monte Carlo method')
pi_Q3 = ((1-0))/numsim * 4 * sum(pi_QMC_Q3)

%% Questtion 4: Newton-Cotes approach: Pythagorean
% Again I use the midpoint rule to compute the integration

% define the width of partition h
h = (1-0)/numsim;

% Define vector of x and y (later filled)
x = zeros(numsim,1);

% x_j = 0 + (j-1/2)h for j=1,...,numsim
for ind = 1:numsim
    x(ind,1) = 0 + (ind- 1/2)*h;
end

% Compute the sequence of values of the function pi_root and approximate
% the integration
display('Problem 3: Pythagorean fomula with Newton-Cotes method')
pi_NC_Q4 =  4 * h * sum(pi_root(x))


%% question 5:

numsim_list= [100,1000,10000];
realpi = pi;

% Dart-Throwing
% Implement numerical integration using QMC with different number of draws
% We simulate 200 times and compute the squared error
DT_ErrPMC_200 = ones(200,3);
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    seed = 1534561;
    for sim = 1:200
        
        seed = seed + sim ;
        rng(seed);
        xy = rand(numsim,2); 
        x = xy(:,1);
        y = xy(:,2);
        pi_DT_QMC = pi_ind(x,y);
        % comute squared residual for QMC
        DT_ErrPMC_200(sim,i) = (realpi - ((1-0)/numsim * 4 * sum(pi_DT_QMC)))^2;
    end
    clear x y
end
% Mean squaed error is obtained as
DT_MErrPMC = sum(DT_ErrPMC_200)/200;

% Pythagorean
% Implement numerical integration using QMC with different number of draws
% We simulate 200 times and compute the squared error
Py_ErrPMC_200 = ones(200,3);
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    seed = 1534561;
    for sim = 1:200
        seed = seed + sim ;
        rng(seed);
        x = rand(numsim,1); 
        % comute squared residual for QMC
        pi_QMC = pi_root(x);
        Py_ErrPMC_200(sim,i) = (realpi - ((1-0))/numsim * 4 * sum(pi_QMC))^2;
    end
    clear x
end
% Mean squaed error is obtained as
Py_MErrPMC = sum(Py_ErrPMC_200)/200;
        

% Dart Throwing
% We then use Newton-Coates
DT_ErrNC = ones(1,3);
% Implement numerical integration using NC and compute the squared error
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    % define the width of partition h
    h = (1-0)/numsim;
    % Define vector of x and y (later filled)
    x = zeros(numsim,1);
    y = zeros(numsim,1);
    
    for ind = 1:numsim
        x(ind,1) = 0 + (ind- 1/2)*h;
        y(ind,1) = 0 + (ind- 1/2)*h;
    end
    pi_NC = ones(numsim,1);
    for ind = 1:numsim
        x_1 = x(ind,1)*ones(numsim,1);
        pi_NC_x = pi_ind(x_1,y);
        pi_NC(ind,1) =  h * sum(pi_NC_x);
    end
    DT_ErrNC(1,i) = (realpi - (4 * h * sum(pi_NC)))^2;
   clear x y  
end
% Mean squaed error is obtained as

% Pythagorean
% We then use Newton-Coates
Py_ErrNC = ones(1,3);
% Implement numerical integration using NC and compute the squared error
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    h = (1-0)/numsim;
    x = zeros(numsim,1);
    for ind = 1:numsim
        x(ind,1) = 0 + (ind- 1/2)*h;
    end
    NC = 4* h * sum(pi_root(x));
    Py_ErrNC(1,i) = (realpi - (4* h * sum(pi_root(x))))^2;
    clear x
end
% Mean squaed error is obtained as


display('Problem 5: Comparison')
display('Dart-Throwing with Monte-Carlo: MSEs (100, 1000, and 10000 draws)')
DT_MErrPMC
display('Pythagorean with Monte-Carlo: MSEs (100, 10000, and 10000 draws)')
Py_MErrPMC

display('Dart-Throwing with Newton-Cotes: squared error (100, 1000, and 10000 nodes)')
DT_ErrNC
display('Pythagorean with Newton-Cotes: squared error (100, 1000, and 10000 nodes)')
Py_ErrNC


diary off