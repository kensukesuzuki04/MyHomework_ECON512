% Spring 2019
% ECON512 Empirical Method
% Homework 7 --- Dynamic Game
% Kensuke Suzuki
% kxs974@psu.edu
% all is good, check plus

delete HW7log.txt
diary('HW7log.txt')
diary on

clear;
tic;
setupParams;

%% Problem 1

% Initial guess
cmat = kron(c, ones(1,L));
p1 = (repmat(c,1,L) + v)./2;
V1 = (p1 - cmat)./(1-beta);

diff = 1;
iter=1;

% solution.mat contains the policy and value function which is solved
% through the algorithm below
% in order to reduce time of computation for each run, 
% I store the result and call it (you may comment out three lines below)
load solution.mat   
p1 = solution.price;
V1 = solution.value;

while diff > 1e-5 && iter < 100000;
    
    % get W
    W = getW(V1);

    % solve for p using fsolves
    psol = @(p) focp(p,W);
    p1_new = fsolve(psol, p1);
    
    % get value function
    V1_new = getV(p1,p1_new,W);
    
    % compute difference
    diff_p = max( abs(((p1_new - p1)./(1+abs(p1_new)))) );
    diff_V = max( abs(((V1_new - V1)./(1+abs(V1_new)))) );
    diff = max([diff_p,diff_V])
    
    % update value and policy functions
    V1 = lambda .* V1_new + (1-lambda).* V1;
    p1 = lambda .* p1_new + (1-lambda).* p1;
    iter = iter +1
    
end

 solution.price = p1;
 solution.value = V1;
 save solution;
 

figure(1);
mesh(V1);
title('Value Function');
saveas(gcf,'value.png')

figure(2);
mesh(p1);
title('Price');
saveas(gcf,'price.png')

%% Problem 2
% transition matrix
% Raw: (omega1, omega2)
% Column (omega1', omega2')
Trans = zeros(L*L,L*L);

for i = 1:L*L
    for j = 1:L*L
        if rem(i,30) == 0
            i1 = fix(i/30);
        else
            i1 = fix(i/30) +1 ;
        end
        i2 = rem(i,30);
        if i1 == 31
            i1 = 30;
        end
        if i2 == 0
            i2 = 30;
        end
        
        if rem(j,30) == 0
            j1 = fix(j/30);
        else
            j1 = fix(j/30) +1 ;
        end
                
        j2 = rem(j,30);
        if j1 == 31
            j1 = 30;
        end
        if j2 == 0
            j2 = 30;
        end
        [i,j];
        ipair = [i1,i2]; % today's state for player 1 and 2
        jpair = [j1,j2]; % future state
        
        
        Trans(i,j)  = Pr(i1,j1,1)*Pr(i2,j2,1)*(1 - D(p1(i1,i2),p1(i2,i1))- D(p1(i2,i1), p1(i1,i2)) ) ...
                    + Pr(i1,j1,2)*Pr(i2,j2,1)*(D(p1(i1,i2),p1(i2,i1)) ) ...
                    + Pr(i1,j1,1)*Pr(i2,j2,2)*(D(p1(i2,i1), p1(i1,i2)) ) ;
                
    end
end



clear state state_new


% initial state
state_int = zeros(1,L*L);
state_int(1,1) = 1;

% 10 period
state = state_int;
for t = 2:10
    state_new =   state * Trans;
    state = state_new;
end

Dstrbn10 = zeros(L,L);
for i = 1:L
    Dstrbn10(i,:) = (state_new((i-1)*L+1:i*L))';
end

% 20 period
state = state_int;
for t = 2:20
    state_new = state * Trans;
    state = state_new;
end

Dstrbn20 = zeros(L,L);
for i = 1:L
    Dstrbn20(i,:) = (state_new((i-1)*L+1:i*L))';
end

% 30 period
state = state_int;
for t = 2:30
    state_new = state * Trans;
    state = state_new;
end

Dstrbn30 = zeros(L,L);
for i = 1:L
    Dstrbn30(i,:) = (state_new((i-1)*L+1:i*L))';
end

figure(3);
mesh(Dstrbn10);
title('Distribution after 10 periods');
saveas(gcf,'10period.png')

figure(4);
mesh(Dstrbn20);
title('Distribution after 20 periods');
saveas(gcf,'20period.png')

figure(5);
mesh(Dstrbn30);
title('Distribution after 30 periods');
saveas(gcf,'30period.png')

%% Problem 3
% stationary distribution

state = state_int;
diff = 1;
while diff > 1e-6
    state_new = state * Trans;
    diff = max(abs(state - state_new))
    state = state_new;
end

StDstrbn = zeros(L,L);
for i = 1:L
    StDstrbn(i,:) = (state_new((i-1)*L+1:i*L))';
end

figure(6);
mesh(StDstrbn);
title('Stationary Distribution');
saveas(gcf,'stationary.png')

diary off