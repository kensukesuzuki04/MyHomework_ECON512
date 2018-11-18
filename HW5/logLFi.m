function fval = logLFi(beta,gamma,u,X,Y,Z)
epsi = (beta * X + gamma * Z + u  ) ;

expepsi = exp(-epsi);
%expepsi = exp(-epsi);
%expepsi = single(expepsi);
logitval =  ( 1 + expepsi ).^(-1);
val1 = Y.* log(logitval) + (ones(20,100)-Y).* log(ones(20,100)-logitval);
val2 = double(val1);
fval = sum(val2);
end
