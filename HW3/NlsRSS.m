function RSS = NlsRSS(beta,X,y)
res = y - exp(X*beta);
RSS = (res' * res);
end
