% ECON512 Homework 3
% Kensuke Suzuki
clear all
delete HW3log.txt
diary('HW3log.txt')
diary on

disp('ECON512 HOMEWORK3: Ken Suzuki')

disp(' ')

data = load('hw3.mat');
X = data.X;
y = data.y;
Parameters = {'beta0'; 'beta1'; 'beta2'; 'beta3';'beta4'; 'beta5'};
%bindcell = cellstr(bind);
%% Problem 1: Nelder Mead Simplex Method

options = optimset('PlotFcns',@optimplotfval, 'Display','iter');

TobitLLF_beta =  @(beta) TobitLLF(beta,X,y);
%beta0 = [1,1,1,1,1,1]';
beta0 = [log(mean(y)),zeros(1,5)]';


%beta_p1 = fminsearch(TobitLLF_beta, beta0, options);
beta_p1 = fminsearch(TobitLLF_beta, beta0);

EstimatedCoeff = beta_p1;
LLF= -TobitLLF(beta_p1,X,y)
disp('------Problem 1------')
disp('Estimated Parameter Vector')
Result_Q1 = table(Parameters, EstimatedCoeff)

%% Problem 2: Quasi-Newton Optimization Method
clear beta 
beta_p2 = fminunc(TobitLLF_beta, beta0)
EstimatedCoeff = beta_p2;

disp('------Problem 2------')
disp('Method: Quasi-Newton Optimization (fminunc)')
disp('Estimated Parameter Vector')
Result_Q2 = table(Parameters, EstimatedCoeff)


%% Problem 3: 

NlsRSS_beta =  @(beta) NlsRSS(beta,X,y);

beta_p3_1 = lsqnonlin(NlsRSS_beta,beta0)
beta_p3_2 = lsqnonlin(NlsRSS_beta,beta_p1)

Coeff1 = beta_p3_1; 
Coeff2 = beta_p3_2;

%% Problem 3: NLS with lsqnonlin

NlsRSS_beta =  @(beta) NlsRSS(beta,X,y);

beta_p3_1 = lsqnonlin(NlsRSS_beta,beta0)
beta_p3_2 = lsqnonlin(NlsRSS_beta,beta_p1)
Coeff1 = beta_p3_1; 
Coeff2 = beta_p3_2;

disp('------Problem 3------')
disp('NLS with lsqnonlin')
disp('Estimated Parameter Vector')
disp('Coeff1: [log(mean(y)),zeros(1,5)] as initial guess')
disp('Coeff2: MLE estimator (question 1) as initial guess')
Result_Q3 = table(Parameters, Coeff1, Coeff2)


%% Problem 4: NLS with Nelder-Mead

beta_p4_1 = fminsearch(NlsRSS_beta,beta0)
beta_p4_2 = fminsearch(NlsRSS_beta,beta_p1)
Coeff1 = beta_p4_1; 
Coeff2 = beta_p4_2;

disp('------Problem 4------')
disp('NLS with Nelder-Mead')
disp('Estimated Parameter Vector')
disp('Coeff1: [log(mean(y)),zeros(1,5)] as initial guess')
disp('Coeff2: MLE estimator (question 1) as initial guess')
Result_Q3 = table(Parameters, Coeff1, Coeff2)

%% Problem 5

intbeta0 = [1:0.1:5];

% MLE with fminsearch
MLE_fmins = zeros(6,size(intbeta0,2));
for i = 1:size(intbeta0,2)
    beta0 = [intbeta0(1,i),zeros(1,5)]';
    beta = fminsearch(TobitLLF_beta, beta0);
    MLE_fmins(:,i)=beta;
end
result_MLE_fmins = [intbeta0; MLE_fmins];

%MLE with quasi newton
MLE_fminunc = zeros(6,size(intbeta0,2));
for i = 1:size(intbeta0,2)
    beta0 = [intbeta0(1,i),zeros(1,5)]';
    beta = fminunc(TobitLLF_beta, beta0);
    MLE_fminunc(:,i)=beta;
end
result_MLE_fminunc = [intbeta0; MLE_fminunc];

%NLS with lsqnonlin
NLS_lsqnonlin = zeros(6,size(intbeta0,2));
for i = 1:size(intbeta0,2)
    beta0 = [intbeta0(1,i),zeros(1,5)]';
    beta = lsqnonlin(NlsRSS_beta,beta0);
    NLS_lsqnonlin(:,i)=beta;
end
result_NLS_lsqnonlin = [intbeta0; NLS_lsqnonlin];

%NLS with lsqnonlin
NLS_fmins = zeros(6,size(intbeta0,2));
for i = 1:size(intbeta0,2)
    beta0 = [intbeta0(1,i),zeros(1,5)]';
    beta = fminsearch(NlsRSS_beta,beta0);
    NLS_fmins(:,i)=beta;
end
result_NLS_fmins = [intbeta0; NLS_fmins];

% MLE with fminsearch
plot(intbeta0, result_MLE_fmins(2,:), ...
    intbeta0, result_MLE_fmins(3,:), ...
    intbeta0, result_MLE_fmins(4,:), ...
    intbeta0, result_MLE_fmins(5,:), ...
    intbeta0, result_MLE_fmins(6,:), ...
    intbeta0, result_MLE_fmins(7,:))
title('Estimated Coefficients: MLE with Nelder-Mead')
xlabel('Initial Guess of beta0')
ylabel('Estimated Coefficients')
legend('beta0','beta1', 'beta3', 'beta5', 'beta5', 'beta6')
saveas(gcf, 'MLE_NM.png')

% MLE with Quasi Newton
plot(intbeta0, result_MLE_fminunc(2,:), ...
    intbeta0, result_MLE_fminunc(3,:), ...
    intbeta0, result_MLE_fminunc(4,:), ...
    intbeta0, result_MLE_fminunc(5,:), ...
    intbeta0, result_MLE_fminunc(6,:), ...
    intbeta0, result_MLE_fminunc(7,:))
title('Estimated Coefficients: MLE with Quasi-Newton')
xlabel('Initial Guess of beta0')
ylabel('Estimated Coefficients')
legend('beta0','beta1', 'beta3', 'beta5', 'beta5', 'beta6')
saveas(gcf, 'MLE_QN.png')

% NLS with lsqnonlin
plot(intbeta0, result_NLS_lsqnonlin(2,:), ...
    intbeta0, result_NLS_lsqnonlin(3,:), ...
    intbeta0, result_NLS_lsqnonlin(4,:), ...
    intbeta0, result_NLS_lsqnonlin(5,:), ...
    intbeta0, result_NLS_lsqnonlin(6,:), ...
    intbeta0, result_NLS_lsqnonlin(7,:))
title('Estimated Coefficients: NLS with lsqnonlin')
xlabel('Initial Guess of beta0')
ylabel('Estimated Coefficients')
legend('beta0','beta1', 'beta3', 'beta5', 'beta5', 'beta6')
saveas(gcf, 'NLS.lsqnonlin.png')

% NLS with fminsearch
figure 
yyaxis left
plot(intbeta0, result_NLS_fmins(2,:))
xlabel('Initial Guess of beta0')
ylabel('Estimated Coefficients')
title('Estimated Coefficients: NLS with Nelder-Mead')
yyaxis right
plot(intbeta0, result_NLS_fmins(3,:), ...
    intbeta0, result_NLS_fmins(4,:), ...
    intbeta0, result_NLS_fmins(5,:), ...
    intbeta0, result_NLS_fmins(6,:), ...
    intbeta0, result_NLS_fmins(7,:))
legend('beta0','beta1', 'beta3', 'beta5', 'beta5', 'beta6')
saveas(gcf, 'NLS_NM.png')

diary off