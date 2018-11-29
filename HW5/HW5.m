% Empirical Method HW5 %
% Ken Suzuki (Penn State)
% kxs974@psu.edu 

clear all
delete HW5log.txt
diary('HW5log.txt')
diary on

% Load Data
load('hw5.mat')

X = data.X;
Y = data.Y;
Z = data.Z;

addpath('../CEtools/');

%% Problem 1
% parameter value
betanot = 0.1;
sigmab = 1;
gamma = 0;

% method
method = 1;

% number of nodes
node = 20;

% set parameter vector
par = [gamma betanot sigmab];

[llf,methodname] = llk_wou(Y,X,Z,par,node,method);

llf = -1 * llf;
% display result
disp('Problem 1')
disp(methodname)
disp('Loglikelihood is')
disp(llf)

%% Problem 2
% parameter value
betanot = 0.1;
sigmab = 1;
gamma = 0;

% method
method = 2;

% number of nodes
node = 100;

% set parameter vector
par = [gamma betanot sigmab];

[llf,methodname] = llk_wou(Y,X,Z,par,node,method);
llf = -1 * llf; % take negative

% display result
disp('Problem 2')
disp(methodname)
disp('Loglikelihood is')
disp(llf)

%% Problem 3

clear par

% number of node
node = 20;

% method: GC
method = 1;

% define function to be minimzied (function of par)
llkwou_min = @(par) llk_wou(Y,X,Z,par,node,method);

% fmincon

% for constraint
A = [0, 0, -1];
b = 0

%[paraGQ, lfGQ] = fminsearch(llkwou_min, [1 1 1] );
[paraGQ, lfGQ] = fmincon(llkwou_min, [1; 1; 1], A, b );
lfGQ = -1 * lfGQ;

% number of node
node = 100;

% method: GC
method = 2;

% define function to be minimzied (function of par)
llkwou_min = @(par) llk_wou(Y,X,Z,par,node,method);

% fminsearch
[paraMC, lfMC] = fmincon(llkwou_min, [1; 1; 1], A, b );
lfMC = -1 * lfMC;


% display result
disp('Problem 3-1 (Gaussian Quadrature)')
disp('Initial guesses are')
disp('   gamma      beta     sigmab')
disp([1 1 1])
disp('   gamma      beta     sigmab')
disp(paraGQ')
disp('Maximized log-likelihood is:')
disp(lfGQ)


% display result
disp('Problem 3-2 (Monte Carlo)')
disp('Initial guesses are')
disp('   gamma      beta     sigmab')
disp([1 1 1])
disp('   gamma      beta     sigmab')
disp(paraMC')
disp('Maximized log-likelihood is:')
disp(lfMC)


%% Problem 4

clear A 
clear b

% method: MC
method = 2;

% number of node
node = 100;

% initial values for parameter vector
gamma = -0.5056;
betanot = 2.5579;
sigmab = 1.1816;
unot = 1;
%sigmaub = 0.9;
rho = 0.9;
sigmau =1;
%intpar = [gamma betanot sigmab unot sigmaub sigmau];
intpar = [gamma betanot sigmab unot rho sigmau];

%define function to be minimized
llkwu_min = @(par) llk_wu(Y,X,Z,par,node,method);

% for constraint
A = [0 0 -1 0 0  0 ; ...
     0 0  0 0 0 -1; ...
     0 0  0 0 -1 0; ...
     0 0  0 0  1 0];
b = [0; 0;  1; 1];

% fmincon
[paraMC, lfMC] = fmincon(llkwu_min, intpar', A, b );

sigmaub = paraMC(3)^(1/2) * paraMC(6)^(1/2) * paraMC(5);
paraMC_cov = paraMC;
paraMC_cov(5) = sigmaub;

sigmaubint = intpar(3)^(1/2) * intpar(6)^(1/2) *intpar(5);
intpar_cov = intpar;
intpar_cov(5) = sigmaubint ;

lfMC = -1 * lfMC;

% display result
disp('Problem 4 (Monte Carlo)')
disp('Initial guesses are')
disp('   gamma      betanot   sigmab    unot      sigmaub   sigmau')
disp(intpar_cov)
disp('Estimated parameters')
disp('   gamma      betanot   sigmab    unot      sigmaub   sigmau')
disp(paraMC_cov')
disp('Maximized log-likelihood is:')
disp(lfMC)

diary off