function fval = getV(p_old,p_new,W) % update new p

global L c beta;

cmat = kron(c, ones(1,L));
D1 = D(p_new, p_old')   ; % firm 1's demand
D2 = D(p_old', p_new)   ; % firm 2's demand 
D0 = ones(L, L) -  D(p_new, p_old')  -  D(p_old', p_new); % demand for the outside good

W0 = W(:,:,1);
W1 = W(:,:,2);
W2 = W(:,:,3);

fval = D(p_new,p_old').*(p_new - cmat) + beta*( D0.* W0 + D1.* W1 + D2.*W2);

end