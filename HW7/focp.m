function fval = focp(p,W) % update new p

global L c beta;

cmat = kron(c, ones(1,L));
D1 = D(p, p')   ; % firm 1's demand
D2 = D(p', p)   ; % firm 2's demand 
D0 = ones(L, L) - D(p, p') - D(p', p); % demand for the outside good

W0 = W(:,:,1);
W1 = W(:,:,2);
W2 = W(:,:,3);

fval = ones(L,L) -  ( ones(1,1) - D1) .* (p - cmat) ... 
        - beta * W1 +  beta * (D0 .* W0 + D1 .* W1 + D2.* W2 );

end