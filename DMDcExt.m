% By Alvin Wang 
% wanghao.bit@gmail.com

function  [Ahat, Bhat, Atilde, Btilde] = DMDcExt(X1, X2, Upsilon, r)
% clasic dmdc when B is unknown
% modified by Alvin Wang from func_DMDc.m by Zhe Bai
    %SVD
    [U,Sig,V] = svd(X1,'econ');

    thresh = 1e-10;
    r = length(find(diag(Sig)>thresh));
    U    = U(:,1:r);
    Sig  = Sig(1:r,1:r);
    V    = V(:,1:r);
    % matrix of the state and input snapshots
    Omega = [X1; Upsilon];
    [Utilde, Stilde, Vtilde] = svd(Omega, 'econ');
    rtilde = length(find(diag(Stilde)>thresh));
    %rtilde =3;
    % truncated r-rank SVD of X2
    [Uhat, ~, ~] = svd(X2, 'econ');

    Utilde1 = Utilde(1:end-rtilde+r, 1:rtilde);
    Utilde2 = Utilde(end-rtilde+r+1:end, 1:rtilde);
    temp = Vtilde(:,1:rtilde)*(Stilde(1:rtilde,1:rtilde))^(-1)*Utilde1';
    Ahat = X2*temp;
    Bhat = X2*Vtilde(:,1:rtilde)*(Stilde(1:rtilde,1:rtilde))^(-1)*Utilde2';
    Atilde = Uhat'*Ahat*Uhat;
    Btilde = Uhat'*Bhat;
end