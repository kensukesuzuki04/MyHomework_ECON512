% Empirical Method HW5 %
% Ken Suzuki (Penn State)
% kxs974@psu.edu 

clear all

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

% fminsearch
[paraGQ, lfGQ] = fminsearch(llkwou_min, [1 1 1] );
lfGQ = -1 * lfGQ;

% number of node
node = 100;

% method: GC
method = 2;

% define function to be minimzied (function of par)
llkwou_min = @(par) llk_wou(Y,X,Z,par,node,method);

% fminsearch
[paraMC, lfMC] = fminsearch(llkwou_min, [1 1 1] );
lfMC = -1 * lfMC;


% display result
disp('Problem 3-1 (Gaussian Quadrature)')
disp('   gamma      beta     sigmab')
disp(paraGQ)
disp('Maximized log-likelihood is:')
disp(lfGQ)


% display result
disp('Problem 3-2 (Monte Carlo)')
disp('   gamma      beta     sigmab')
disp(paraMC)
disp('Maximized log-likelihood is:')
disp(lfMC)


%% Problem 4

% method: MC
method = 2;

% number of node
node = 100;

% initial values for parameter vector
gamma = 1;
betanot = 1;
sigmab = 1;
unot = 1;
sigmaub = 0.5;
sigmau =1;
intpar = [gamma betanot sigmab unot sigmaub sigmau];

%define function to be minimized
llkwu_min = @(par) llk_wu(Y,X,Z,par,node,method);

% fminsearch
[paraMC, lfMC] = fminsearch(llkwu_min, intpar );
lfMC = -1 * lfMC;

% display result
disp('Problem 4 (Monte Carlo)')
disp('   gamma      betanot   sigmab    unot      sigmaub   sigmau')
disp(paraMC)
disp('Maximized log-likelihood is:')
disp(lfMC)
