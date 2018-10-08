function LLF = TobitLLF(beta,X,y)
f = -exp(X*beta) + y .* (X*beta) - log(factorial(y)) ;
%f = -exp(X*beta) + y .* (X*beta) ;
LLF = - sum(f);
end
