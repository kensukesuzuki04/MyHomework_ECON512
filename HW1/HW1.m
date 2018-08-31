
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