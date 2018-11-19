function [llf,methodname] = llk_wou(Y,X,Z,par,node,method)
    % compute negative of llf

gamma = par(1);
betanot = par(2);
sigmab = par(3);
%unot = par(4);
unot = 0;
%sigmaub = par(5);
sigmaub = 0;
%sigmau = par(6);
sigmau = 0;

mu = [betanot unot];
Sigma = [sigmab sigmaub; sigmaub sigmau];

ui = 0;


if method == 1
    methodname = 'Gaussian Quadrature';
        % if method is Gaussian Quadrature
    [rcoef,w] = qnwnorm(node, betanot, sigmab);


for i = 1:length(rcoef)
    betai = rcoef(i,1);
        % pick ith draw of beta
    epsi = (betai * X + gamma * Z + ui  ) ;
    logitval = ( 1 + exp(-1 * epsi) ).^(-1);
        % compute the logistic CDF
    lkt = logitval.^Y .* (ones(20,100)-logitval).^(ones(20,100)-Y);
        % compute the contribution of each year
    lkii(i,:) = prod(lkt);
        % product over years
end
lki = w' * lkii;
    % numerical integration
    
llki = log(lki);

llf = -1 * sum(llki,2);

elseif method == 2
    methodname = 'Monte Carlo';
        % if method is MC
        
    norm = haltonNormShuffle(node, 1, 3);
    rcoef = repmat(betanot,node,1) + sigmab * norm';

for i = 1:length(rcoef)
    betai = rcoef(i,1);
        % pick ith draw of beta
    epsi = (betai * X + gamma * Z + ui  ) ;
    logitval = ( 1 + exp(-1 * epsi) ).^(-1);
        % compute the logistic CDF
    lkt = logitval.^Y .* (ones(20,100)-logitval).^(ones(20,100)-Y);
        % compute the contribution of each year
    lkii(i,:) = prod(lkt);
        % product over years
end
lki = sum(1/node * lkii);
    % numerical integration
    
llki = log(lki);

llf = -1 * sum(llki,2);        

end

    
end
