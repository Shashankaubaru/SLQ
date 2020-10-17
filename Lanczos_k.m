function [T,V,w] = Lanczos_k(A,v0,k)
%
% k-step basic Lanczos iteration
%
% usage:    [T,V,w] = Lanczos_k(A,v0,k)
t  = norm(v0,2);
  v  =  v0/t; 
  beta = 0; 
  V = [v]; 
  vold = v;
%%-------------------- for stopping -- 
  orthTol = 1.e-08;
  wn = 0.0 ;
%%-------------------- main loopLanSVD.m
 for j=1:k
    w = A*v; 
    w = w - beta*vold ; 
    alpha = w'*v; 
    wn = wn + alpha*alpha;
    T(j,j) = alpha; 
    w = w - alpha*v;
%%-------------------- full reorthogonalization
    t1 = V'*w;
%%fprintf(1,' it %4d inPr %e \n',k,norm(t));
    w = w - V* t1;     
    beta = w'*w;
    if (beta*j < orthTol*wn) 
        break
    end           
%%-------------------- for orth. test
    wn   = wn+2.0*beta;
    beta = sqrt(beta) ;
    vold = v; 
    v = w/ beta;
    V(:,j+1) = v;
    T(j,j+1) = beta; 
    T(j+1,j) = beta; 
 end


