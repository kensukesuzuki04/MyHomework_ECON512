function val = pi_ind(x,y)

ind = x.^2 + y.^2;
val = zeros(length(x),1);

for sim = 1:length(x)
if ind(sim,1) <= 1
    val(sim,1) = 1;
else
    val(sim,1) = 0;
end
end

