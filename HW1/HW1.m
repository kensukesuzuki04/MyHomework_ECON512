
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
F = A(2,:)

% Solve linear equatuons
x = inv(A)*b

%% Problem 4

