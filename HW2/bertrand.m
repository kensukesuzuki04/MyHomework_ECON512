function fval = bertrand(p,v)

% for given vector of p and v, solve demand
fval = exp(v - p) ./ ( 1 + sum( exp(v - p) ) );

end