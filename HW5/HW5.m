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
% Number of node
node = 20;

% parameter value
betanot = 0.1;
sigmab = 1;
gamma = 0;

% Taking u out of the model
ui = 0;
% 
% tic
% logLFGQ = -1 * logLFGQ(betanot,gamma,X,Y,Z,ui, sigmab, node)
% toc
% 
% % display result
% disp('Problem 1: Loglikelihood is')
% disp(logLFGQ)

%% Problem 2
% Number of node
node = 100;

% parameter value
betanot = 0.1;
sigmab = 1;
gamma = 0;

% Taking u out of the model
ui = 0;

tic
logLFMC_P2 = -1* logLFMC(betanot,gamma,X,Y,Z,ui,sigmab, node);
toc

% display result
disp('Problem 2: Loglikelihood is')
disp(logLFMC_P2)

%% Problem 3

%logLFGQ_min = @(beta_vect) logLFGQ(beta_vect(1),gamma,X,Y,Z,beta_vect(2),node);

% Use Gaussian 
%x = fminsearch(@(beta_vect) logLFGQ(beta_vect(1),gamma,X,Y,Z,beta_vect(2),node), [0.1 1] )

labP3 = ['beta0' 'sigma_b'  'gamma'];

disp('Problem 3-1: Estimated parameter is:')
disp(labP3)
disp('Maximized log-likelihood is:')
disp(' ')

node = 100;
% use MC
logLFMC_min = @(para) logLFMC(para(1,1), para(1,2), X, Y, Z, ui, para(1,3), node);
[estpara, minLF] = fminsearch(logLFMC_min, [0.1 0 1]);
logLFMC = -1 * minLF;

% display result
disp('Problem 3-2:  Estimated parameter is:')
disp('    beta     gamma     sigmab')
disp(estpara)
disp('Maximized log-likelihood is:')
disp(logLFMC)


%% Problem 4

node = 100;
% use MC
logLFMCu_min = @(para) logLFMCu(para(1,1), para(1,2), X, Y, Z, para(1,3), para(1,4), para(1,5), para(1,6), node);

[estpara, minLFu] = fminsearch(logLFMCu_min, [0.1 0.1 0.1 0.5 0.5 1]);
logLFMCu = -1 * minLFu;