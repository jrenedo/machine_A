%% Calculate the rotor losses:
% 
clc
clear all



% Modified bessel function of the first kind: I = besseli(nu,Z)
% Modified bessel function of the second kind: K = besselk(nu,Z)
% k_p=sqrt(j*omega*mu*sigma)




%% Parameters (Ali machine 1): 

mu_0=4*pi*10^-7; % [m kg s^-2 A^-2]
p=2; % pairs of poles


h_space=5; % space order
q_eval=p*h_space; % space order
k_space=6; % time order

B_given=50*10^-3; % [T]

n_mech=65000; % mechanical speed rpm
f_1=n_mech/60*p; % [Hz]
omega=f_1*2*pi*k_space;




R_1=20*10^-3; % [m] rotor radius
R_2=27.4*10^-3; % [m] magnet radius
R_3=31.15*10^-3; % [m] stator radius
R_4=50*10^-3; % [m] outer radius
t_sl=2.95*10.0^-3; % [m] sleeve thickness
R_wave=R_2+t_sl;

L=109*10^-3; % [m] axial length

sigma_1_eval=6.7*10^6; % [S/m] 
sigma_2_eval=0.77*10^6; % [S/m] 
sigma_3_eval=3*10^-15; % [S/m] 
sigma_4_eval=3*10^-15; % [S/m] 


mu_1=750*mu_0; % [m kg s^-2 A^-2]
mu_2=1.07*mu_0; % [m kg s^-2 A^-2]
mu_3=mu_0; % [m kg s^-2 A^-2]
mu_4=5000*mu_0; % [m kg s^-2 A^-2]

%% Define the symbolic variable

% syms r
% q=1;
% 
% y=besseli(q,r)

syms sigma_1
syms sigma_2
syms sigma_3
syms sigma_4

% sigma_1=6.7*10^6; % [S/m] 
% sigma_2=0.77*10^6; % [S/m] 
% sigma_3=3*10^-15; % [S/m] 
% sigma_4=3*10^-15; % [S/m] 




%% Current sheet model:


syms q

J_kq=1;
B1=[0 0 0 0 0 J_kq*mu_4];
B=transpose(B1)


k_p_1=sqrt(j*omega*mu_1*sigma_1)
k_p_2=sqrt(j*omega*mu_2*sigma_2)
k_p_3=sqrt(j*omega*mu_3*sigma_3)
k_p_4=sqrt(j*omega*mu_4*sigma_4)

a_11=besseli(q,k_p_1*R_1)
a_12=-besseli(q,k_p_2*R_1)
a_13=-besselk(q,k_p_2*R_1)

a_21=mu_2/mu_1*k_p_1*besseli_d(q,k_p_1*R_1)
a_22=-k_p_2*besseli_d(q,k_p_2*R_1)
a_23=-k_p_2*besselk_d(q,k_p_2*R_1)

a_32=besseli(q,k_p_2*R_2)
a_33=besselk(q,k_p_2*R_2)
a_34=-besseli(q,k_p_3*R_2)
a_35=-besselk(q,k_p_3*R_2)

a_42=-mu_3/mu_2*k_p_2*besseli_d(q,k_p_2*R_2)
a_43=-mu_3/mu_2*k_p_2*besselk_d(q,k_p_2*R_2)
a_44=k_p_3*besseli_d(q,k_p_3*R_2)
a_45=k_p_3*besselk_d(q,k_p_3*R_2)

a_54=-besseli(q,k_p_3*R_3)
a_55=-besselk(q,k_p_3*R_3)
a_56=besseli(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk(q,k_p_4*R_3);


a_64=-mu_4/mu_3*k_p_3*besseli_d(q,k_p_3*R_3)
a_65=-mu_4/mu_3*k_p_3*besselk_d(q,k_p_3*R_3)
a_66=k_p_4*(besseli_d(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk_d(q,k_p_4*R_3))




syms A

A(1,1)=a_11;
A(1,2)=a_12;
A(1,3)=a_13;

A(2,1)=a_21;
A(2,2)=a_22;
A(2,3)=a_23;

A(3,2)=a_32;
A(3,3)=a_33;
A(3,4)=a_34;
A(3,5)=a_35;

A(4,2)=a_42;
A(4,3)=a_43;
A(4,4)=a_44;
A(4,5)=a_45;

A(5,4)=a_54;
A(5,5)=a_55;
A(5,6)=a_56;


A(6,4)=a_64;
A(6,5)=a_65;
A(6,6)=a_66;




%A2=transpose(A)

%% Estimation of the coefficients


%C=inv(A);

X_1=A\B;

C1_temp=X_1(1);
C_1=double(subs(C1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

C2_temp=X_1(2);
C_2=double(subs(C2_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

D2_temp=X_1(3);
D_2=double(subs(D2_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

C3_temp=X_1(4);
C_3=double(subs(C3_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

D3_temp=X_1(5);
D_3=double(subs(D3_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

C4_temp=X_1(6);
C_4=double(subs(C4_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

D4_temp=-X_1(6)*besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4);
D_4=double(subs(D4_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));


% Until here it seems to be working.

%% Calculate the power losses for J_kq=1:

S1=2*pi*R_1*L;
S2=2*pi*R_2*L;

P1_j1_temp=1/2*real(-j*omega*C_1*besseli(q,k_p_1*R_1)/mu_1*conj(-k_p_1*C_1*besseli_d(q,k_p_1*R_1)))*S1;
P1_j1=double(subs(P1_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

P2_j1_temp=1/2*real(-j*omega*(C_2*besseli(q,k_p_2*R_2)+D_2*besselk(q,k_p_2*R_2))/mu_2*conj(-k_p_2*C_2*besseli_d(q,k_p_2*R_2)+D_2*besselk_d(q,k_p_2*R_2)))*S2;
P2_j1=double(subs(P2_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));

P_j1_mag=P2_j1-P1_j1
P_j1_hub=P1_j1

%% Solution with no eddy currents:

syms q



clear A
% q=q_eval;
% A=zeros(7);
syms A


A(1,1)=R_1^q;
A(1,2)=-R_1^q;
A(1,3)=-R_1^(-q);

A(2,2)=R_2^q;
A(2,3)=R_2^(-q);
A(2,4)=-R_2^q;
A(2,5)=-R_2^(-q);

A(3,4)=R_3^q;
A(3,5)=R_3^(-q);
A(3,6)=-R_3^q;
A(3,7)=-R_3^(-q);

A(4,6)=R_4^q;
A(4,7)=R_4^(-q);

A(5,1)=-1/mu_1*R_1^(q);
A(5,2)=1/mu_2*R_1^(q);
A(5,3)=-1/mu_2*R_1^(-q);

A(6,2)=-1/mu_2*R_2^(q);
A(6,3)=1/mu_2*R_2^(-q);
A(6,4)=1/mu_3*R_2^(q);
A(6,5)=-1/mu_3*R_2^(-q);

A(7,4)=-1/mu_3*R_3^(q);
A(7,5)=1/mu_3*R_3^(-q);
A(7,6)=1/mu_4*R_3^(q);
A(7,7)=-1/mu_4*R_3^(-q);

B1=[0 0 0 0 0 0 J_kq/q_eval*R_3];
B=transpose(B1);

X_1=A\B;

C10_temp=X_1(1);
C_10=double(subs(C10_temp,q, q_eval))

C20_temp=X_1(2);
C_20=double(subs(C20_temp,q, q_eval))

D20_temp=X_1(3);
D_20=double(subs(D20_temp,q, q_eval))

C30_temp=X_1(4);
C_30=double(subs(C30_temp,q, q_eval))

D30_temp=X_1(5);
D_30=double(subs(D30_temp,q, q_eval))

C40_temp=X_1(6);
C_40=double(subs(C40_temp,q, q_eval))

D40_temp=X_1(7);
D_40=double(subs(D40_temp,q, q_eval))


B_calculated=-q_eval*(C_20*R_wave^(q_eval-1)+D_20*R_wave^(-q_eval-1))

K_B=B_given/B_calculated

%% Final estimation of the rotor losses:

P_mag=P_j1_mag*K_B^2
P_hub=P_j1_hub*K_B^2

P_tot=P_mag+P_hub

