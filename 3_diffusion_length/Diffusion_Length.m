clear all; clc; close all


%% Step 1.
%beta_ and tau_s (ns)
beta_s=  0.6192;  
tau_s=168.4;
%% Step 2.
%diffusion coefficient (nm2/ns)
D=978.3639

%% extract tau_e (ns) interpolating time at 1/e intensity
time=(0:0.1:250);
tau_e=interp1(exp(-(time./tau_s).^beta_s),time,1/exp(1))

%% diffusion length (nm)
diff_length=sqrt(tau_e*D)
