function [llf,methodname] = llk_wu(Y,X,Z,par,node,method)
    % compute negative of llf

gamma = par(1);
betanot = par(2);
sigmab = par(3);
unot = par(4);
sigmaub = par(5);
sigmau = par(6);


mu = [betanot unot];
Sigma = [sigmab sigmaub; sigmaub sigmau];
U = chol(Sigma);

if method == 2
    methodname = 'Monte Carlo';
        % if method is MC
        
    norm = haltonNormShuffle(node, 2, 2);
    rcoef = repmat(mu,node, 1) + (U' * norm)';

for i = 1:length(rcoef)
    betai = rcoef(i,1);
    ui = rcoef(i,2);
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
