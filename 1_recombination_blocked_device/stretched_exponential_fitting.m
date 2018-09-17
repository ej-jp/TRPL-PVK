clear all; clc; close all

%load experimental data 
load block_device_data.csv;% time (ns) vs normalized counts (-)
time=block_device_data(:,1);
counts=block_device_data(:,2);




fit_stretchexp = fit( time, counts, 'exp(-(x./tau_s)^beta_s)', 'Startpoint', [0.9,202])
figure(1)
semilogy(time,feval(fit_stretchexp,time), time, counts,'r+' )
%result
%      General model:
%      fit_stretchexp(x) = exp(-(x./tau_s)^beta_s)
%      Coefficients (with 95% confidence bounds):
%        beta_s =      0.6192  (0.6006, 0.6379)
%        tau_s =       168.4  (164.9, 171.8)

% uncomment to save data for plots
%data=[time,feval(fit_stretchexp,time), time, counts];
%csvwrite('data_fit.out',data);

