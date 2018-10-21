% Empirical method HW4
% Kensuke Suzuki
% Penn State
% October 20

clear all
delete HW4log.txt
diary('HW4log.txt')
diary on

disp('ECON512 HOMEWORK4: Ken Suzuki')

%% Questtion 1: Quasi-Monte Carlo method

% Number of random draw
numsim = 10000;

seed = 1534561;
rng(seed);

% Define function: pi_ind
% returns 1 if x^2 + y^2 <= 1 and 0 otherwise

% Generate random sequence for x and y using rand
seq = rand(numsim,2);
x = seq(:,1);
y = seq(:,2);

% We now compute the sequence of values of indicator function using the
% random sequence generated above
pi_QMC_Q1 = pi_ind(x,y);

% Compute numerical integation 
display('Problem 1: Quasi-Monte Carlo method')
pi_Q1 = ((1-0)*(1-0))/numsim * 4 * sum(pi_QMC_Q1)

clear x y
%% Questtion 2: Newton-Cotes approach
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


%% Questtion 3: Newton-Cotes approach: another functional form
% I use Halton sequence to generate random draws

% Define function: pi_root
% returns (1-x^2)^(1/2)

seed = 1534561;
rng(seed);
% Generate random sequence for x 
x = rand(numsim,1);

% We now compute the sequence of values of indicator function using the
% random sequence generated above
pi_QMC_Q3 = pi_root(x);

% Compute numerical integation 
display('Problem 3: Pythagorean fomula with Quasi-Monte Carlo method')
pi_Q3 = ((1-0))/numsim * 4 * sum(pi_QMC_Q3)

%% Questtion 4: Newton-Cotes approach: another functional form
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

numsim_list= [1000,10000,100000];
realpi = pi;



% Implement numerical integration using QMC with different number of draws
% We simulate 200 times and compute the squared error
ErrQMC_200 = ones(200,3);
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    seed = 1534561;
    for sim = 1:200
        seed = seed + sim ;
        rng(seed);
        % comute squared residual for QMC
        x = rand(numsim,1);
        pi_QMC = pi_root(x);
        ErrQMC_200(sim,i) = (realpi - ((1-0))/numsim * 4 * sum(pi_QMC))^2;
    end
    clear x
end
% Mean squaed error is obtained as
MErrQMC = sum(ErrQMC_200)./numsim_list;
        

% We then use Newton-Coates
ErrNC = ones(1,3);
% Implement numerical integration using NC and compute the squared error
for i = 1:length(numsim_list)
    numsim = numsim_list(1,i);
    h = (1-0)/numsim;
    x = zeros(numsim,1);
    for ind = 1:numsim
        x(ind,1) = 0 + (ind- 1/2)*h;
    end
    NC = 4* h * sum(pi_root(x));
    ErrNC(1,i) = (realpi - (4* h * sum(pi_root(x))))^2;
    clear x
    
end
% Mean squaed error is obtained as
ErrNC;

display('Problem 5: Comparison')
display('Monte-Carlo mean squared error (1000, 10000, and 100000 draws)')
MErrQMC
display('Newton-Cotes squared error (1000, 10000, and 100000 nodes)')
ErrNC

diary off