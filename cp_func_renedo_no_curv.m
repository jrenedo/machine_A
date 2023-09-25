function [ output_args ] = cp_func_renedo( Rm,Rs,Rr,bo,Qs,R,dslot )
%Calculates the CP according to Renedo.
%   Detailed explanation goes here


%----------------------------------------------------------------------
%----------------------------------------------------------------------
mu0=4*pi*1e-7;
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%-------------------------------------------------------------------------
Mtype='P'; %Type of magnetization
Nlambda=64; %Number of Fourier coefficients in the relative air gap
%permeance distribution (max. 128)
NB=301; %Number of Fourier coefficients in the flux density distribution
N=300; %number of evaluation points (always use an even number!)
%----------------------------------------------------------------------
%Polar coordinates of the characteristic points on the slot outline
Rp1=Rr;
Rp2=Rs;
Rp3=Rs;
Rp4=Rs;
thetas=2*pi/Qs; %slot pitch in rad
alphao=2*asin(0.5*bo/Rs); %slot opening in rad
theta1=thetas/2-alphao/2;
theta2=theta1+alphao;
%----------------------------------------------------------------------

%---------------------

% Try to reproduce eaxctly the same as ali's code but with the SC toolbox
% for the derivative. 

path(path,'C:\local\Documents Jaime\PhD\MATLAB\Driscoll SA toolbox\sc')


Rg=(Rr+Rs)/2;

N_points=1000;

% Parameters of the machine:

% alpha_t=tau_w/Rs;
% alpha_s=2*pi/Qs-alpha_t;

alpha_s=bo/Rs;
alpha_t=2*pi/Qs-alpha_s;


t=Rg*alpha_t;
s=Rg*alpha_s;
d=4*dslot;
delta=R-Rr;
g2=Rs-Rr;

theta_lambda=alpha_t+alpha_s;
theta_points=(0:(N_points))*theta_lambda/N_points;

% Coordinates of the problem:


z_dom= [Rr Rs Rs+j*t/2 Rs+d+j*t/2 Rs+d+j*(t/2+s) Rs+j*(t/2+s) Rs+j*(t+s) Rr+j*(t+s)];
z_points=Rr+delta+j*(0:N_points)*(t+s)/N_points;





p=polygon(z_dom);
% Indicates the right angles in the Canonical Domain:
alpha=[0.5 0.5 1 1 1 1 0.5 0.5];

% Remember that acording to this criteria X is the plane with the toothed
% member and w is the plane with the canonical rectangle.

% Define the Canonical Domain:

f=crrectmap(p,alpha);



% Vertices of the canonical rectangle:
vc1=evalinv(f,z_dom(1));
vc2=evalinv(f,z_dom(2));
vc7=evalinv(f,z_dom(7));
vc8=evalinv(f,z_dom(8));
% figure
% plot([vc1 vc2 vc7 vc8])

g_poly=vc1-vc2;

l_poly=abs(vc1-vc8);
% Modified rectangle:


tp_points=evalinv(f,z_points);
%SC_derivative=evaldiff(f,evalinv(f,z_points(2:(N_points-1))))
for k=2:(N_points)
    %SC_derivative(k)=evaldiff(f,evalinv(f,z_points(k)));
    SC_derivative(k)=(evaldiff(f,tp_points(k)))^-1;
end
SC_derivative(1)=(evaldiff(f,evalinv(f,z_points(2))))^-1;
SC_derivative(N_points+1)=(evaldiff(f,evalinv(f,z_points(N_points-1))))^-1;

%d_tz=real(SC_derivative)*g2/g_poly+j*imag(SC_derivative)*(s+t)/l_poly;
d_tz=real(SC_derivative)*g2/g_poly+j*imag(SC_derivative)*g2/g_poly;
%SC_derivative2=evaldiff(f,tp_points(2:(N_points-1)))

%t_points=real(tp_points)*log(Rr/Rs)/g_poly+j*imag(tp_points)*theta_lambda/l_poly+log(Rr);
% plot(t_points,'r')

%k_points=exp(t_points);

for k=1:(N_points+1)
    d_ks(k)=d_tz(k);
end



cp_func_nc=-conj(d_ks);



% OLD 



output_args=[theta_points; cp_func_nc];


end

