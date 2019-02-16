function fval = getW(V1) % update new p

global L Pr;

W = zeros(L,L,3);
for i = 1:1:L
    for j=1:1:L
        W(i,j,1) = Pr(i,:,1) * (V1 * Pr(j,:,1)'); % W0
        W(i,j,2) = Pr(i,:,2) * (V1 * Pr(j,:,1)'); % W1
        W(i,j,3) = Pr(i,:,1) * (V1 * Pr(j,:,2)'); % W2
    end
end

fval = W;

end