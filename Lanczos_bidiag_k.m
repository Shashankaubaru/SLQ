function [B,U,V] = Lanczos_bidiag_k(A,v0,k)
%
% k-step Golub-Kahan-Lanczos bidiagonalization iteration
%
% usage:   [T,V,w] = Lanczos_bidiag_k(A,v0,k)
t  = norm(v0,2);
v  =  v0/t;
beta = 0;
V = [v];
uold =zeros(size(A,1),1);
%%-------------------- for stopping --
orthTol = 1.e-08;
wn = 0.0 ;
%%-------------------- main loopLanSVD.m
for j=1:k
    u = A*v;
    u = u - beta*uold ;
    %%-------------------- full reorthogonalization of U
    if(j>1)
        t2=U'*u;
        u = u - U* t2;
    end
    alpha = norm(u);
    u=u/alpha;
    U(:,j) = u;
    wn = wn + alpha*alpha;
    B(j,j) = alpha;
    v = A'*u - alpha*v;
    %%-------------------- full reorthogonalization of V
    t1 = V'*v;
    v = v - V* t1;
    beta = v'*v;
    if (beta*j < orthTol*wn)
        break
    end
    %%-------------------- for orth. test
    wn   = wn+2.0*beta;
    beta = sqrt(beta) ;
    uold = u;
    v = v/ beta;
    if(j<k)
        V(:,j+1) = v;
        B(j,j+1) = beta;
    end
end


