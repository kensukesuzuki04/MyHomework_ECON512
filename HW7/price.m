function fval = price(p,W) % update new p

global L c beta ;

cmat = kron(c, ones(1,L));
D1 = D(p, p')   ; % firm 1's demand
D2 = D(p', p)   ; % firm 2's demand 
D0 = ones(L, L) - D(p, p') - D(p', p); % demand for the outside good

W1 = W(:,:,1);
W2 = W(:,:,2);
W3 = W(:,:,3);

fval = cmat + (1 - beta * W2 + beta * (D0 .* W1 + D1 .* W2 + D2.* W3 )) ./ ( 1 - D1);


end