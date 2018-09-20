% ECON512 Homework 2 
% Kensuke Suzuki
clear all
delete HW2log.txt
diary('HW2log.txt')
diary on

disp('ECON512 HOMEWORK2: Ken Suzuki')

%% Define bertrand and bertrandfoc function
% bertrand: return demand for each good
% bertrandfoc: return system of FOC (LHS)

%% Problem 1

p = [1;1];
v = [2;2];

Ans1 = bertrand(p,v)

D0 = 1 / (1+ sum(exp(v-p)) );

P1 = sprintf('Problem1: for vA=vB=2 and pA=pB=1, DA= %f, DA= %f, and D0= %f.', Ans1(1,1),Ans1(2,1),D0);
disp(P1);

%% Problem 2

clear all

v = [2;2];

p = [1;1];

fVal_foc = @(p) bertrandfoc(p,v);
i_fVal_foc = fVal_foc(p);

iJac = eye(size(p,1));

maxit = 100;
tol = 1e-6;

tic
for iter = 1:maxit
    fnorm = norm(i_fVal_foc);
    fprintf('iter %d: p(1)=%f, p(2)=%f, norm(f(x))=%.8f \n', iter, p(1),p(2),norm(i_fVal_foc));
    if norm(i_fVal_foc)<tol
        break
    end
    d = -(iJac*i_fVal_foc);
    p = p + d;
    fOld_foc = i_fVal_foc;
    i_fVal_foc = fVal_foc(p);
    u = iJac*(i_fVal_foc - fOld_foc);
    iJac = iJac + ( (d-u)* (d'*iJac) )/(d'*u);
end
elapsedTime_p2 = toc;

P2 = sprintf('Problem2: for vA=vB=2, equilibrium prices are: PA= %f, PB= %f; time elapsed is %f.', p(1,1),p(2,1), elapsedTime_p2);
disp(P2);



%% Problem 3

clear all

v = [2;2];
p = [1;1];

fVal_foc = @(p) bertrandfoc(p,v);
fVal_focg = @(p,g) bertrandfocg(p,v,g);

maxit = 100;
tol = 1e-6;

tic
for iter = 1:maxit
    
    fval = fVal_foc(p);
    if norm(fval) < tol
        break
    end
    fprintf('iter %d: p(1)=%f, p(2)=%f, norm(f(x))=%.8f \n', iter, p(1),p(2),norm(fval));
    
    % set pOld
    pOld = [2;2];
    pA_Old = pOld(1,1); 
    
    % compute the LHS of FOC for good A price
    fOld_1 = fVal_focg(pOld,1);    
    
    % for given pB, solve the first equation for pA
    % We use Secant Method
    for iter_1 =1:maxit
        fval_1 = fVal_focg(p,1);
        if abs(fval_1) < tol
            break
        else
            pA_New = p(1,1) - ( (p(1,1) - pA_Old) / (fval_1 - fOld_1) )* fval_1;
            pA_Old = p(1,1);
            p(1,1) = pA_New;
            fOld_1 = fval_1;
        end
    end   
    
    % Use the solution for pA obtained above, solve for pB
    pB_Old = pOld(2,1); 
    fOld_2 = fVal_focg(pOld,2);    
    for iter_2 = 1:maxit
        fval_2 = fVal_focg(p,2);
        if abs(fval_2) < tol
            break
        else
            pB_New = p(2,1) - ( (p(2,1) - pB_Old) / (fval_2 - fOld_2) )* fval_2;
            pB_Old = p(2,1);
            p(2,1) = pB_New;
            fOld_2 = fval_2;
        end
    end
    
end
elapsedTime_p3 = toc;

P3 = sprintf('Problem3: for vA=vB=2, equilibrium prices are: PA= %f, PB= %f; time elapsed is %f.', p(1,1),p(2,1), elapsedTime_p3);
disp(P3);

%% Problem 4

clear all

v = [2;2];

p = [1;1];

fVal_bertrand = @(p) bertrand(p,v);
fVal_foc = @(p) bertrandfoc(p,v);
i_fVal_foc = fVal_foc(p);

maxit = 100;
tol = 1e-6;

tic
for iter = 1:maxit
    fnorm = norm(i_fVal_foc);
    fprintf('iter %d: p(1)=%f, p(2)=%f, norm(f(x))=%.8f \n', iter, p(1),p(2),norm(i_fVal_foc));
    if norm(i_fVal_foc)<tol
        break
    end
    p_next = 1./ ([1;1] - fVal_bertrand(p) );
    p = p_next;
    i_fVal_foc = fVal_foc(p);
end
elapsedTime_p4 = toc;

P4 = sprintf('Problem4: for vA=vB=2, equilibrium prices are: PA= %f, PB= %f; time elapsed is %f.', p(1,1),p(2,1), elapsedTime_p4);
disp(P4);

%% Problem 5

clear all

vB_5 = [0:.2:3];
v_5 = [2*ones(1,size(vB_5,2));vB_5 ];
result = [vB_5; ones(1,size(vB_5,2)); ones(1,size(vB_5,2)) ];

for vindex = 1:size(vB_5,2)
    
    p = [1;1];
    v = v_5(:,vindex);
    
    fVal_foc = @(p) bertrandfoc(p,v);
    i_fVal_foc = fVal_foc(p);
    iJac = eye(size(p,1));
    
    maxit = 100;
    tol = 1e-6;
    
    for iter = 1:maxit
        fnorm = norm(i_fVal_foc);
        %fprintf('iter %d: p(1)=%f, p(2)=%f, norm(f(x))=%.8f \n', iter, p(1),p(2),norm(i_fVal_foc));
        if norm(i_fVal_foc)<tol
            break
        end
        d = -(iJac*i_fVal_foc);
        p = p + d;
        fOld_foc = i_fVal_foc;
        i_fVal_foc = fVal_foc(p);
        u = iJac*(i_fVal_foc - fOld_foc);
        iJac = iJac + ( (d-u)* (d'*iJac) )/(d'*u);
    end
    result(2,vindex) = p(1);
    result(3,vindex) = p(2);
end

plot(vB_5,result(2,:),vB_5,result(3,:))
title('Equilibrium prices along with change in V_B')
xlabel('v_B')
ylabel('Equilibrium Price: P_A,P_B')
legend('P_A','P_B')

diary off