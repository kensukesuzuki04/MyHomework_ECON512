function fval = bertrandfoc(p,v)

% for given vector of p and v, solve demand
D = exp(v - p) ./ ( 1 + sum( exp(v - p) ) );

% First order condition boils down to
fval = ones(size(p,1),1) - diag(ones(size(p,1),1)-D)*p;

end