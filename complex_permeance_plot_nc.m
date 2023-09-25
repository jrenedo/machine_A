% Plot complex permeance function:
% Ali Machine 2:

clc
clear all


%% Parameters:

mu_0=4*pi*10^-7; % [m kg s^-2 A^-2]

R_r=40/2*10^-3; % [m] rotor radius
R_s=62.3/2*10^-3; % [m] stator radius
R_m=54.8/2*10^-3; % [m] magnet radius
bo=3e-3; % [m] slot opening width

p=2; % pole pairs
alpha_p=1; % pole-arc to pitch ratio
Qs=12; % number of slots

c_gap=R_s-R_m; % [m] clearance gap
dslot=10e-3; %slot depth
 

B_r=1.05; % [T] magnet remanence CHECK THIS WITH SULEIMAN
mu_r=1.07; % relative permeability
H_c=B_r/mu_0; % A/m


N_harm=1551; % maximum harmonic order, has to be even
N_points=10000; % number of points for the waveform

r=R_m+.5*10^-3; % [m] radius where we want to calculate the waveform (used in 2DFEA)
r=27*10^-3; % [m] radius where we want to calculate the waveform (used in 2DFEA)
R_wave=r;




N_pos=555;


%% Complex permeance:

% curvature vs no curvature:
%load('cp_function_slotting')


theta=(0:N_points)*2*pi/p/(N_points); % electrical angle


temp2=cp_func_renedo( R_m,R_s,R_r,bo,Qs,R_wave,dslot );
cp_func=temp2(2,:);
theta_points=temp2(1,:); % mechanical angle






%% no curvature

% Calculate the coeffs of the fourier series single slots:
clear f1_t f1_r

temp=cp_func_renedo_no_curv( R_m,R_s,R_r,bo,Qs,R_wave,dslot );
lambda_nc=temp(2,:);
angle_nc=temp(1,:); % electrical degrees
%% 

cp_func2=[cp_func cp_func cp_func];
lambda_nc2=[lambda_nc lambda_nc lambda_nc];

theta_p2=[theta_points-2*pi/Qs theta_points theta_points+2*pi/Qs];
angle_nc2=[angle_nc-2*pi/Qs angle_nc angle_nc+2*pi/Qs];

figure
subplot(2,1,1)
hold on
plot(theta_points/pi*p,real(cp_func),'LineWidth',1.8,'Color',[0    0.4470    0.7410])
plot(angle_nc/pi*p,real(lambda_nc),'-.','LineWidth',1.8,'Color',[0.6350    0.0780    0.1840])
ylabel('Real part, Re$ \{ \lambda  \}$','Interpreter','LaTex')
grid
subplot(2,1,2)
hold on
plot(theta_points/pi*p,-imag(cp_func),'LineWidth',1.8,'Color',[0    0.4470    0.7410])
plot(angle_nc/pi*p,-imag(lambda_nc),'-.','LineWidth',1.8,'Color',[0.6350    0.0780    0.1840])
grid
ylabel('Imaginary part, Im$ \{ \lambda  \}$','Interpreter','LaTex')



xlabel('$\times$ $\pi$, electrical angle $\theta$ (rad)','Interpreter','LaTex')

figure
subplot(2,1,1)
hold on
plot(theta_p2/pi*p,real(cp_func2),'LineWidth',1.8,'Color',[0    0.4470    0.7410])
plot(angle_nc2/pi*p,real(lambda_nc2),'-.','LineWidth',1.8,'Color',[0.6350    0.0780    0.1840])
ylabel('Real part, Re$ \{ \lambda  \}$','Interpreter','LaTex')
grid
xlim([theta_p2(1) theta_p2(end)]/pi*p)
ylim([0.97 1.01])

subplot(2,1,2)
hold on
plot(theta_p2/pi*p,imag(cp_func2),'LineWidth',1.8,'Color',[0    0.4470    0.7410])
plot(angle_nc2/pi*p,imag(lambda_nc2),'-.','LineWidth',1.8,'Color',[0.6350    0.0780    0.1840])
grid
ylabel('Imaginary part, Im$ \{ \lambda  \}$','Interpreter','LaTex')
xlim([theta_p2(1) theta_p2(end)]/pi*p)
legend('CP with curvature','CP no curvature','northeast')

xlabel('$\times$ $\pi$, electrical angle $\theta$ (rad)','Interpreter','LaTex')


%Now if you want a 3 by 2 inc figure:
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 7.5 8])

%To see the result on screen:
set(gcf,'Units','centimeters','Position',[2 2 7.5 8])
