function fval = logLFMCu(betanot, gamma, X, Y, Z, unot, sigmab, sigmaub, sigmau, node)
seed = 3123;
rng(seed);
x = mvnrnd([betanot unot], [sigmab sigmaub; sigmaub sigmau], node);
[sigmab sigmaub; sigmaub sigmau]
for i = 1:length(x)
    betai = x(i,1);
    ui = x(i,2);
    val1(i,:) =  logLFi(betai,gamma,ui,X,Y,Z);
end
val2 = ones(1,node)*(1/node) * val1;
fval = -1 *  sum(val2, 2);


end
