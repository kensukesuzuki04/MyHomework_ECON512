function fval = logLFMC(betanot,gamma,X,Y,Z,ui,sigmab, node)
seed = 3123;
rng(seed);
x = normrnd(betanot,sigmab, [node,1]);

for i = 1:length(x)
    betai = x(i,1);
    val1(i,:) =  logLFi(betai,gamma,ui,X,Y,Z);
end
val2 = ones(1,node)*(1/node) * val1;
fval = -1 *  sum(val2, 2);


end
