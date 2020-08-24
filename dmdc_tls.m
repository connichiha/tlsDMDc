% By Alvin Wang 
% wanghao.bit@gmail.com

function [ Ahat, Bhat ] = dmdc_tls(Data,Upsilon,r)
%%%%%%%
% Script to implement tls dynamic mode decomposition w Control
% This implementation assumes a time series of data with uniform dt
%
% Modified by Alvin Wang from dmd_tls by Scott Dawson
%
%%%%%%%
H1 = Data(:,1:end-1);
H2 = Data(:,2:end);
Omega = [H1;Upsilon];
%% TLS method
[U,S,V] = svd([Omega', H2'],'econ');
V12 = V(1:r,(r+1):end);
V22 = V((r+1):end,(r+1):end);

if rank(V22)<2
    error('TLS solution does not exist')
end

Admdc = -(V12*V22^(-1))'; 
Ahat = Admdc(1:(r-1),1:(r-1));
Bhat = Admdc(1:(r-1),r:end);

%% method in paper: Dawson, 2016, same as above
% [U,S,V] = svd([H1;Upsilon;H2],'econ');
% U11 = U(1:r,1:r);
% U21 = U((r+1):end,1:r);
% Admdc1  = (U21*U11^-1);
% Ahat = Admdc1(:,1:r);
% Bhat = Admdc1(:,r+1:end);
%%
% [Evecs, DMD_eigs] = eig(Admd);
% DMD_modes = Ur*Evecs;
end