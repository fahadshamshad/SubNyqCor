% Reconstract matrix m of rank r and size nXn with SVLS algorithm.
% Input:
% Br - matrix of row measurements (Br = Ar*m + noise*N(0,1))
% Bc - matrix of column measurements (Bc = m*Ac + noise*N(0,1))
% Ar - row design matrix
% Ac - column design matrix
% r - rank of unknown matrix (optional) 
%
% Output:
% X - estimate of M where we minimize  ||m*ac-Bc||^2 +  ||ar*M-Br||^2
% est_r - estimated rank (when rank r is unknown)
%
function [X, est_r] = SVLS( Br,Bc,Ar,Ac,r)

est_r = r;

%%%% Taking SVD %%%%%%

[uC s1 v1] = svd(Bc);
%[u2 s2 uR] = svd(Br);

%%%% truncating u2 and v1  %%%%

uC = uC(:,1:r);
%uR = uR(:,1:r);

% [uC, ~, ~] = lansvd(Bc,r,'L','OPTIONS'); %find the r largest vectors and singulr values of Bc.
% [uR, ~, ~] = lansvd(Br',r,'L','OPTIONS'); %find the r largest vectors and singulr values of Br.
%  Xr=(uR*(Ac'*uR\Bc'))';% find Xr
 Xc=uC*(Ar*uC\Br); % find Xc
%Xr = Bc*pinv(uR'*Ac)*uR';
%Xc = uC*pinv(Ar*uC)*Br;

X = Xc;
%[u ,d, v]  = lansvd(X,r,'L','OPTIONS');
%X = u(:,1:r)*d(1:r,1:r)*v(:,1:r)';
end

