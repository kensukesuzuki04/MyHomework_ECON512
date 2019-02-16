function fval = D(p1,p2)

global v
fval = exp(v-p1) ./ (1 + exp(v-p1) + exp(v-p2));

end