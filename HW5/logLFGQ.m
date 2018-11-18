function fval = logLFGQ(betanot,gamma,X,Y,Z,ui, sigmab, node)
[x,w] = qnwnorm(node, betanot, sigmab);
u=0;

for i = 1:length(x)
    betai = x(i,1);
    val1(i,:) =  logLFivpa(betai,gamma,ui,X,Y,Z);
end
val2 = w' * val1;
fval = -1 * sum(val2, 2);

end
