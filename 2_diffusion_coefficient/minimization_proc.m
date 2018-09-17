 
clear all; clc; close all

%% Parameters previously obtained in blocked device (step 1)
%beta_ and tau_s (ns)
global beta_s 
global tau_s
beta_s=0.6192;  
tau_s=168.4;

%% minimization procedure
% coeffiecient D=2000 nm^2/ns is equivalent to 0.02 cm2/s
% theta is the vector of optimized parameters. theta0 is the initial value
% and thetamax, thetamin limit the range of permitted values. 
thetamax=[1200];
theta0=  [978.3639];
thetamin=[800];



%[D,fval,exitflag,output]=fminsearchcon(@main_diffusion_funk,theta0,thetamin,thetamax,[],[],[],optimset('Algorithm','Interior-point','TolFun',1e-2,'TolX',1e-2,'MaxFunEvals',100,'MaxIter',100))


[D]=fminsearchcon(@main_diffusion_funk,theta0,thetamin,thetamax,[],[],[],optimset('Algorithm','Interior-point','TolFun',1e-2,'TolX',1e-2,'MaxFunEvals',100,'MaxIter',100))


%confidence interval for D coefficient
% bootstraping.m