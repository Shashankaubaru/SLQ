function [spnorm,running_avg, sum_vl] = Lanc_Quad_Schatten_norm( X, m, nvecs, p)
%% function ld = Lanc_Quad_Schatten_norm( A, m, nvecs)
% The function computes an approximate Schatten p-norm of any 
% (rectangluar or square) matrix X  using the Stochastic Lanczos Quadrature

%-- Inputs
% X - the  input matrix
% m - Number of Lanczos steps (degree)
% nvecs - Number of starting vectors
% p - value of p for the Schatten p-norm norm

%-- Output
% spnorm - The Schatten p-norm of X estimated by SLQ
% running_avg - Running average of the estimated norm 
% sum_vl - Individual estimates for each starting vector v_l


%%- By Shashanka Ubaru
%% Initialization

 n = size(X,2);
cnt_est=0;
running_avg = zeros(nvecs,1);
sum_vl =  zeros(nvecs,1);
%% Main loop
 for ii = 1:nvecs
   w = sign(randn(n,1)); % Random radamacher vector
    v0 = w /norm(w);
    [B,~,~] = Lanczos_bidiag_k(X,v0,m);
    [~,D,V1]=svd(B);
    theta  = abs(diag(D));
    gamma2 = V1(1,:).^2;

    %% sum of gamma2*theta^p
    thetap = theta.^p;
    count=sum(gamma2'.*(thetap));
    sum_vl(ii) = (count)*n;
    cnt_est = (cnt_est+count);
    running_avg(ii) = n*(cnt_est/ii);
 end
spnorm=running_avg(end);