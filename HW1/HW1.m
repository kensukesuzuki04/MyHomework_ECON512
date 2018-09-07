% ECON512 Homework 1 
% Kensuke Suzuki
clear all
delete HW1log.txt
diary('HW1log.txt')
diary on

disp('ECON512 HOMEWORK1: Ken Suzuki')
disp(' ')

%% Problem 1
X = [1,1.5,3,4,5,7,9,10];
Y1 = -2 + .5*X;
Y2 = -2 + 0.5 * X.^2;

plot(X, Y1, X, Y2)
title('Problem 1: plot Y1 and Y2 against X')
xlabel('X') % x-axis label
ylabel('Y1 and Y2') % y-axis label
legend('Y1 = -2 + .5X','Y2 = -2 + .5X^2')

disp('------Problem 1------')
X
Y1
Y2
disp('Figure is saved in the directory')
disp(' ')

%% Problem 2
% Create 200x1 vector X 
clear X
X = linspace(-10,20,200)';
sumX = sum(X);

disp('------Problem 2------')
disp('First define the vector X using linspace. (display of X is ommitted)')
disp('Sum of all elements in X is')
sumX
disp(' ')


%% Problem 3
A = [2,4,6; 1,7,5; 3,12,4];
b = [-2;3;10];

% C
C = A'*b;
% D
D = inv(A'*A) * b;

% E
E0 = A .* (b*ones(1,3));
E = sum(sum(E0),2);

% F
F0 = [A(1,:);A(3,:)];
F = [F0(:,1), F0(:,2)];

% Solve linear equatuons
x = inv(A)*b;
% it's better to stick with the A\b even though new matlab will do exactly
% the same thing


disp('------Problem 3------')
A
b
C
D
E
F
x
disp(' ')

%% Problem 4
% block diagonal matrix
B =  blkdiag(A,A,A,A,A);
% use kron() instead

disp('------Problem 4------')
B
disp(' ')

%% Problem 5
clear A

rng(2);
A = normrnd(10,5,[5,3]);

disp('------Problem 5------')
disp('Before converting with 1 and 0 matrices')
A

for i = 1:size(A,1)
    for j = 1:size(A,2)
        if A(i,j) < 10
            A(i,j) = 0;
        else
            A(i,j) = 1;
        end
    end
end

disp('After converting with 1 and 0 matrices')
A
disp(' ')

%% Problem 6
clear X
clear Y

filename = 'datahw1.csv';
data = csvread(filename);
% csvread replaces NaN with 0. it creates bias in the estimates

X = [ones(4392,1), data(:,3), data(:,4), data(:,6)];
Y = data(:,5);

% Pointe estimates
betahat = inv(X'*X)*X'*Y;

% Standard error
% residual
e = Y - (X * betahat);
sigmahat = (e'* e)/(size(X,1)-size(X,2));
cov = sigmahat * inv(X'*X);
var = diag(cov);
stderr = var.^(1/2);

betase = [betahat';stderr'];
rowNames = {'Coefficient','Std Error'};
colNames = {'beta0', 'beta1', 'beta2', 'beta3'};


disp('------Problem 6------')
disp('Estimated Coefficients and Standard Errors')
Table = array2table(betase,'RowNames',rowNames,'VariableNames',colNames);
disp(Table)

disp('end of HW1')

%tab = [betaind; tab0]

diary off