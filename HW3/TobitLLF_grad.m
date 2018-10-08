function [LLF, grad] = TobitLLF_grad(beta,X,y)
f = -exp(X*beta) + y .* (X*beta) - log(factorial(y)) ;
grad = -(-X' * exp(X*beta) +  X' * y);
%f = -exp(X*beta) + y .* (X*beta) ;
LLF = - sum(f);
end
