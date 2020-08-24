% https://github.com/eurika-kaiser/SINDY-MPC/blob/master/utils/
% By Alvin Wang 
% wanghao.bit@gmail.com

%function [sysmodel_DMDc,U,Up] = DMDc_SS(X,U,dt)
function [Ar, Br, U,Up] = DMDc_SS(X,U)
X1 = X(:,1:end-1);
X2 = X(:,2:end);
numVar = size(X1,1); numOutputs = numVar; numInputs = size(U,1);
Gamma = U;%U(1:end-1);
Omega = [X1; Gamma];
[U,S,V] = svd(Omega,'econ'); 
G = X2*V*S^(-1)*U';
[Up,Sp,Vp] = svd(X2,'econ'); 
Ar = G(:,1:numVar);
Br = G(:,numVar+1:end);
%Cr = eye(numVar); 
%sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);