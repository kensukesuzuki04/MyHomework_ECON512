
% ECON512 Homework 1 
% Kensuke Suzuki

clear all
%% Problem 1

X = [1,1.5,3,4,5,7,9,10];
Y1 = -2 + .5*X;
Y2 = -2 + 0.5 * X.^2;

plot(X, Y1, X, Y2)
title('Problem 1: plot Y1 and Y2 against X')
xlabel('X') % x-axis label
ylabel('Y1 and Y2') % y-axis label
legend('Y1 = -2 + .5X','Y2 = -2 + .5X^2')

%% Problem 2
% Create 200x1 vector X 
clear X
X = linspace(-10,20,200)';
sumX = sum(X)

%% Problem 3

A = [2,4,6; 1,7,5; 3,12,4]
b = [-2;3;10]

% C
C = A'*b
% D
D = inv(A'*A) * b

% E
E0 = A .* (b*ones(1,3));
E = sum(sum(E0),2)

% F
F0 = [A(1,:);A(3,:)]
F = [F0(:,1), F0(:,2)]

% Solve linear equatuons
x = inv(A)*b

%% Problem 4

% block diagonal matrix
B =  blkdiag(A,A,A,A,A);

%% Problem 5
clear A

A = normrnd(10,5,[5,3])

for i = 1:size(A,1)
    for j = 1:size(A,2)
        if A(i,j) < 10
            A(i,j) = 0;
        else
            A(i,j) = 1;
        end
    end
end

disp(A)

%% Problem 6

filename = 'datahw1.csv';
data = csvread(filename);
indp = [ones(4392,1), data(:,3), data(:,4), data(:,6)];
dpn = data(:,5);

% Pointe estimates
betahat = inv(indp'*indp)*indp'*dpn

% Standard error
% residual
e = dpn - (indp * betahat);
sigmahat = (e'* e)/(size(indp,1)-size(indp,2));
cov = sigmahat * inv(indp'*indp);
var = diag(cov);
stderr = var.^(1/2)



