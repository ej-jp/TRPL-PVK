% main routine created by ej-jp 
% Date 9 Jan 2015
% 1D diffusion model 
% Version 1.1
% Model Description
% This model describes the dynamics of the decaying population of carriers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SSE = mainx8(theta)

D=theta(1);


global beta_s 
global tau_s

%% Initialization
N=500; % number of points to discretize the layer thickness
L=491;  % total length of the perovskite layer, nm

% absortion coefficient of perovskite nm^-1 (cm^-1 divided by 1e7)
% this value is obtained from absorption measurement
alpha=2.10e4/1e7;
n0=1;%arbitrary intensity at x=0 and t=0


% x-vector containing N points
x=(0:(L-0)/(N-1):L);

%illumination side comment and uncomment as requerided
% illumination from left side
n=n0*exp(-alpha*(x));
%lightning from right side
%n=n0*exp(-alpha*(L-x));

t0 = 0.0001;%ns
tf = 250;%ns
tspan = [t0 tf];
init0 =n;


% The boundaries for the pvk layer are algrbraic equations and 
% are not treated as differential equations. The way to
% let ode15s differ between the two kinds of equations is by
% introducing the mass matrix, M:
M=eye(N,N);
% First and last point (boundaries) in the vector are algebraic
% equations
M(1,1)=0.0;
M(N,N)=0.0;
%The matrix is set as sparse
M=sparse(M);
%set options
options = odeset('Mass',M);
%% 

function deriv=diffusion(t,n_vari)
% all x nodes
ni=n_vari(1:N);
%inner x nodes
ni_i=ni(2:N-1);

%expression for k(t) as derivative of the stretched exponential
k=beta_s*tau_s^(-beta_s)*(t).^(beta_s-1);


%% Discretation, finds the partial derivatives in x-direction
%Estimates the first derivatives of ni
dnidx = dss020(0,L,N,ni,1)';
%Estimates the second derivatives of ni
d2nidx2 = dss042(0,L,N,ni,dnidx,1,2)';
% Diffusion equation model itself
dnidt=D*d2nidx2(2:N-1)-k.*ni_i;

% Boundary conditions.
% Dirichlet boundary condition (perfect quenching) at x=0
boundary_1 = ni(1);
% Neumann boundary condition (perfect blocked side) at x=L
boundary_2 = dnidx(N);
% Merging all the vectors into one column vector
deriv=[boundary_1;dnidt;boundary_2];
end

[t,nn] = ode15s(@diffusion,tspan,init0,options);


nxt=nn(:,1:N);


%%%SSE evaluation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration of population across the layer thickness for each time
for i=1:size(t)
it(i)=trapz(x,nxt(i,1:N));
end
it_norm=it/max(it);

%load experimental data 
load quenched_device_data.csv;% time (ns) vs normalized counts (-)
time=quenched_device_data(:,1);
counts=quenched_device_data(:,2);



%Spline functions of experimental and theoretical data for easy
%substration
t_eval=(t0:0.1:tf);
spline_exp1=spline(time,counts,t_eval);
spline_mod1=spline(t,it_norm,t_eval);

% evaluation of SSE
% sum of squared errors
 SSE=sum((spline_exp1-spline_mod1)).^2;
% alternatively, squared sum of logarithm difference
% SSE=sum((log(spline_exp1)-log(spline_mod1))).^2;


% live view of the optimization process
figure(1)
hh=semilogy(time,counts,'rx',t,it_norm);
set(hh(2),'linewidth',3);



% post-procession upon D estimation. Hold commented during optimization
%mesh plot of n(x,t)
% figure(3)
% mesh(x,t,nxt)
% 
% %export data 
% a=t;
% b=it_norm';
% c=time;
% d=counts;
% data_mod=[a,b];
% data_exp=[c,d];
% csvwrite('data_mod.out',data_mod);
% csvwrite('data_exp.out',data_exp);
















end

