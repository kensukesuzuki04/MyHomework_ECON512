clear;
global L l kappa rho eta c ...
    delta beta  CRIT lambda v Omega Pr;

L = 30;
l = 15;
kappa = 10;
rho = 0.85;

eta = log(0.85)/log(2);

c = zeros(L,1);
c(1:l) = kappa .* [1:l]'.^eta;
c(l+1:end) = kappa * l^eta;

delta = 0.03;

v = 10;
beta = 1/1.05;

CRIT = 1e-10; %10^(-6);
lambda = 0.9; %.75; %For Dampening

Omega = (1:1:30)';
%Omega = omega';

% Transition matrix
Pr = zeros(L,L,2); % 
for h =1:2 % q=0 for dimension 1
    q = h-1;
    for i = 1:1:L
        for j = 1:1:L
            if j == i + q
                Pr(i,j,h) = 1 - (1-(1-delta)^i);
            elseif j == i + q -1
                Pr(i,j,h) = 1-(1-delta)^i;
            else
                Pr(i,j,h) = 0;
            end
        end
    end
end

Pr(1,1,1) = 1;
Pr(1,2,1) = 0;
Pr(L,L,2) = 1;