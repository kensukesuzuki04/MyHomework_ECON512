% Spring 2019
% ECON512 Empirical Method
% Homework 6 --- Dynamic Programing
% Kensuke Suzuki
% kxs974@psu.edu

clear all;

%% Basic setup
delta = 0.95;
p0 = 0.5;
rho = 0.5;
sigma = 0.1;

N = 201 ;% number of grid
k = linspace(0,100,N);

%% Problem 1
disp('Refer to the PDF')

%% Problem 2 --- Tauchen

Z = 21; % grid points generated

[prob,grid]=tauchen(Z,p0,rho,sigma);
disp(['The dimensions of prob are ' num2str(size(prob)) ])
disp(['The dimensions of grid are ' num2str(size(grid)) ])


%% Problem 3 --- Value function iteration

vinitial=zeros(N,Z);
vrevised=zeros(N,Z);
decision=zeros(N,Z);

invest=kron(ones(1,Z),k');
disp(['The dimensions of harvest  are ' num2str(size(invest)) ])
ONE=ones(N,1);

diff = 1;
while diff > 1e-6
    
    Ev = vinitial*prob' ;% (NxZ * ZxZ = NxZ);
    
    Ev;
    
    for i = 1:N
        harvest = ( kron(ones(N,Z),k(:,i)') - invest) ;
        harvest_z = max(harvest, 0);
        pi = (ones(N,1)*grid) .* harvest - 0.2 * (harvest_z.^(1.5));
        check = pi + delta*Ev;
        [vrevised(i,:),decision(i,:)] = max( pi + delta*Ev );
        decision;
    end
    diff = norm(vrevised - vinitial) / norm(vrevised);
    vinitial = vrevised;
    
end

%% Problem 4

derule=zeros(N,Z);

for i=1:Z
    derule(:,i)=k(decision(:,i))';
end


plot(k,derule)
hold on
plot(k,k)
xlabel('current capital');
ylabel('next period capital');
title('optimal investment decision')

%% Problem 5 Expected path

T = 1:20;

path_k = zeros(3,length(T));
path_p = zeros(3,length(T));
ind = find(grid==1)
path_k(:,1) = [NaN; N; NaN]; % Nth column corresponds k = 100
path_p(:,1) = [NaN; ind; NaN]; % 11th column corresponds to p = 1

for t = 2:20
    path_k(2,t) = decision(path_k(2,t-1), decision(path_p(2,t-1)));
    exp_p = prob(path_p(2,t-1),:)*grid';
    diff_p = 1./abs(grid-exp_p);
    [a,b]=  max(diff_p);
    path_p(2,t) = b;
end



%% Problem 6 Tauchen with coarse grid

clear prob grid

Z = 5; % grid points generated

[prob,grid]=tauchen(Z,p0,rho,sigma);
disp(['The dimensions of prob are ' num2str(size(prob)) ])
disp(['The dimensions of grid are ' num2str(size(grid)) ])

vinitial=zeros(N,Z);
vrevised=zeros(N,Z);
decision=zeros(N,Z);

invest=kron(ones(1,Z),k');
disp(['The dimensions of harvest  are ' num2str(size(invest)) ])
ONE=ones(N,1);

diff = 1;
while diff > 1e-6
    
    Ev = vinitial*prob' ;% (NxZ * ZxZ = NxZ);
    
    Ev;
    
    for i = 1:N
        harvest = ( kron(ones(N,Z),k(:,i)') - invest) ;
        harvest_z = max(harvest, 0);
        pi = (ones(N,1)*grid) .* harvest - 0.2 * (harvest_z.^(1.5));
        check = pi + delta*Ev;
        [vrevised(i,:),decision(i,:)] = max( pi + delta*Ev );
        decision;
    end
    diff = norm(vrevised - vinitial) / norm(vrevised);
    vinitial = vrevised;
    
end


derule=zeros(N,Z);

for i=1:Z
    derule(:,i)=k(decision(:,i))';
end

figure
plot(k,derule)
hold on
plot(k,k)
xlabel('current capital');
ylabel('next period capital');
title('optimal investment decision')
