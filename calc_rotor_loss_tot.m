%% Calculate the rotor losses:
%
clc
clear all
close all


% Modified bessel function of the first kind: I = besseli(nu,Z)
% Modified bessel function of the second kind: K = besselk(nu,Z)
% k_p=sqrt(j*omega*mu*sigma)
%



%% Parameters (Padova):

mu_0=4*pi*10^-7; % [m kg s^-2 A^-2]
p=2; % pairs of poles


n_rpm=65000; % [rpm]

f_1=n_rpm/60*p; % [Hz]


alpha_p=1;

R_1=40/2*10^-3; % [m] rotor radius
R_2=54.8/2*10^-3; % [m] magnet radius
R_3=62.3/2*10^-3; % [m] stator radius
R_4=108/2*10^-3; % [m] outer radius

L=109*10^-3; % [m] axial length

sigma_1_eval=6.7*10^6; % [S/m] 
sigma_2_eval=0.77*10^6; % [S/m] 
sigma_3_eval=3*10^-15; % [S/m] 
sigma_4_eval=3*10^-15; % [S/m] 



mu_1=750*mu_0; % [m kg s^-2 A^-2]
mu_2=1.07*mu_0; % [m kg s^-2 A^-2]
mu_3=mu_0; % [m kg s^-2 A^-2]
mu_4=5000*mu_0; % [m kg s^-2 A^-2]

% delta_1=sqrt(2/(sigma_1_eval*omega*mu_1)); % Skin depth in [m]
% delta_2=sqrt(2/(sigma_2_eval*omega*mu_2)); % Skin depth in [m]
% delta_3=sqrt(2/(sigma_3_eval*omega*mu_3)); % Skin depth in [m]
% delta_4=sqrt(2/(sigma_4_eval*omega*mu_4)); % Skin depth in [m]



load('amplitude_harmonics')



%% Loop to calculate the total losses:



%% Select only the significant harmonics:

B_treshold=0.001;
vec_space2=[];
vec_time2=[];
vec_B_given=[];

B_treshold_n=0.0001;
vec_space2_n=[];
vec_time2_n=[];
vec_B_given_n=[];

% Mat_multiple_slot_export=Mat_multiple_slot_export*0;
% Mat_multiple_slot_export(4,2)=1;

for count1=2:length(vec_time)
    for count2=1:length(vec_space)
        if Mat_multiple_slot_export(count2,count1)>=B_treshold
            vec_B_given=[vec_B_given Mat_multiple_slot_export(count2,count1)];
            vec_time2=[vec_time2 vec_time(count1)]
            vec_space2=[vec_space2 vec_space(count2)]
        end
        omega=2*pi*f_1*vec_time(count1);
        vec_s_depth1(count1)=sqrt(2/(sigma_1_eval*omega*mu_1))*1000; % [mm]
        vec_s_depth2(count1)=sqrt(2/(sigma_2_eval*omega*mu_2))*1000; % [mm]

    end
    
    if Mat_multiple_slot_n_export(count2,count1)>=B_treshold
            vec_B_given_n=[vec_B_given_n Mat_multiple_slot_n_export(count2,count1)];
            vec_time2_n=[vec_time2_n vec_time(count1)];
            vec_space2_n=[vec_space2_n vec_space(count2)];
            
        end

end

%% Calculation of rotor losses forward rotating:

length_harmonics=length(vec_B_given);

P_ms_matrix=zeros(length_harmonics,1);

P_ms_tot=0;
    for count=1:length_harmonics
        k_time=vec_time2(count)
        
        h_space=vec_space2(count);
        q_eval=p*h_space; % space order
        omega=f_1*2*pi*k_time;
        B_ms_given=vec_B_given(count)
        
        
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
        B=transpose(B1);
        
        
        k_p_1=sqrt(j*omega*mu_1*sigma_1);
        k_p_2=sqrt(j*omega*mu_2*sigma_2);
        k_p_3=sqrt(j*omega*mu_3*sigma_3);
        k_p_4=sqrt(j*omega*mu_4*sigma_4);
        
        a_11=besseli(q,k_p_1*R_1);
        a_12=-besseli(q,k_p_2*R_1);
        a_13=-besselk(q,k_p_2*R_1);
        
        a_21=mu_2/mu_1*k_p_1*besseli_d(q,k_p_1*R_1);
        a_22=-k_p_2*besseli_d(q,k_p_2*R_1);
        a_23=-k_p_2*besselk_d(q,k_p_2*R_1);
        
        a_32=besseli(q,k_p_2*R_2);
        a_33=besselk(q,k_p_2*R_2);
        a_34=-besseli(q,k_p_3*R_2);
        a_35=-besselk(q,k_p_3*R_2);
        
        a_42=-mu_3/mu_2*k_p_2*besseli_d(q,k_p_2*R_2);
        a_43=-mu_3/mu_2*k_p_2*besselk_d(q,k_p_2*R_2);
        a_44=k_p_3*besseli_d(q,k_p_3*R_2);
        a_45=k_p_3*besselk_d(q,k_p_3*R_2);
        
        a_54=-besseli(q,k_p_3*R_3);
        a_55=-besselk(q,k_p_3*R_3);
        a_56=besseli(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk(q,k_p_4*R_3);
        
        
        a_64=-mu_4/mu_3*k_p_3*besseli_d(q,k_p_3*R_3);
        a_65=-mu_4/mu_3*k_p_3*besselk_d(q,k_p_3*R_3);
        a_66=k_p_4*(besseli_d(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk_d(q,k_p_4*R_3));
        
        
        
        
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
        
        P1_j1_temp=1/2*real(-1i*omega*C_1*besseli(q,k_p_1*R_1)/mu_1*conj(-k_p_1*C_1*besseli_d(q,k_p_1*R_1)))*S1;
        P1_j1=double(subs(P1_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));
        
        P2_j1_temp=1/2*real(-1i*omega*(C_2*besseli(q,k_p_2*R_2)+D_2*besselk(q,k_p_2*R_2))/mu_2*conj(-k_p_2*(C_2*besseli_d(q,k_p_2*R_2)+D_2*besselk_d(q,k_p_2*R_2))))*S2;
        P2_j1=double(subs(P2_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));
        
        P_j1_mag=P2_j1-P1_j1;
        P_j1_hub=P1_j1;
        
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
        C_10=double(subs(C10_temp,q, q_eval));
        
        C20_temp=X_1(2);
        C_20=double(subs(C20_temp,q, q_eval));
        
        D20_temp=X_1(3);
        D_20=double(subs(D20_temp,q, q_eval));
        
        C30_temp=X_1(4);
        C_30=double(subs(C30_temp,q, q_eval));
        
        D30_temp=X_1(5);
        D_30=double(subs(D30_temp,q, q_eval));
        
        C40_temp=X_1(6);
        C_40=double(subs(C40_temp,q, q_eval));
        
        D40_temp=X_1(7);
        D_40=double(subs(D40_temp,q, q_eval));
        
        
%         B_calculated=-q_eval*(C_20*R_2^(q_eval-1)+D_20*R_2^(-q_eval-1));
        B_calculated=abs(-j*q_eval*(C_30*R_wave^(q_eval-1)+D_30*R_wave^(-q_eval-1)));
        
        
        % CHANGES HERE
        K_B_ms=B_ms_given/B_calculated;
        
        
        %% Final estimation of the rotor losses:
        
        % Multiple slots:
        P_mag_ms=alpha_p*P_j1_mag*K_B_ms^2;
        P_hub_ms=alpha_p*P_j1_hub*K_B_ms^2;
        
        P_ms_matrix(count)=P_mag_ms;
        P_ms_tot=P_ms_tot+P_mag_ms+P_hub_ms;
        
        
        
        
    end


%% Calculation of rotor losses backward rotating:

length_harmonics=length(vec_B_given_n);

P_ms_matrix_n=zeros(length_harmonics,1);

P_ms_tot_n=0;
    for count=1:length_harmonics
        k_time=vec_time2_n(count);
        
        h_space=vec_space2_n(count);
        q_eval=p*h_space; % space order
        omega=f_1*2*pi*k_time;
        B_ms_given=vec_B_given_n(count);
        
        
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
        B=transpose(B1);
        
        
        k_p_1=sqrt(j*omega*mu_1*sigma_1);
        k_p_2=sqrt(j*omega*mu_2*sigma_2);
        k_p_3=sqrt(j*omega*mu_3*sigma_3);
        k_p_4=sqrt(j*omega*mu_4*sigma_4);
        
        a_11=besseli(q,k_p_1*R_1);
        a_12=-besseli(q,k_p_2*R_1);
        a_13=-besselk(q,k_p_2*R_1);
        
        a_21=mu_2/mu_1*k_p_1*besseli_d(q,k_p_1*R_1);
        a_22=-k_p_2*besseli_d(q,k_p_2*R_1);
        a_23=-k_p_2*besselk_d(q,k_p_2*R_1);
        
        a_32=besseli(q,k_p_2*R_2);
        a_33=besselk(q,k_p_2*R_2);
        a_34=-besseli(q,k_p_3*R_2);
        a_35=-besselk(q,k_p_3*R_2);
        
        a_42=-mu_3/mu_2*k_p_2*besseli_d(q,k_p_2*R_2);
        a_43=-mu_3/mu_2*k_p_2*besselk_d(q,k_p_2*R_2);
        a_44=k_p_3*besseli_d(q,k_p_3*R_2);
        a_45=k_p_3*besselk_d(q,k_p_3*R_2);
        
        a_54=-besseli(q,k_p_3*R_3);
        a_55=-besselk(q,k_p_3*R_3);
        a_56=besseli(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk(q,k_p_4*R_3);
        
        
        a_64=-mu_4/mu_3*k_p_3*besseli_d(q,k_p_3*R_3);
        a_65=-mu_4/mu_3*k_p_3*besselk_d(q,k_p_3*R_3);
        a_66=k_p_4*(besseli_d(q,k_p_4*R_3)-besseli(q,k_p_4*R_4)/besselk(q,k_p_4*R_4)*besselk_d(q,k_p_4*R_3));
        
        
        
        
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
        
        P1_j1_temp=1/2*real(-1i*omega*C_1*besseli(q,k_p_1*R_1)/mu_1*conj(-k_p_1*C_1*besseli_d(q,k_p_1*R_1)))*S1;
        P1_j1=double(subs(P1_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));
        
        P2_j1_temp=1/2*real(-1i*omega*(C_2*besseli(q,k_p_2*R_2)+D_2*besselk(q,k_p_2*R_2))/mu_2*conj(-k_p_2*C_2*besseli_d(q,k_p_2*R_2)+D_2*besselk_d(q,k_p_2*R_2)))*S2;
        P2_j1=double(subs(P2_j1_temp,[q, sigma_1, sigma_2, sigma_3, sigma_4], [q_eval, sigma_1_eval, sigma_2_eval, sigma_3_eval, sigma_4_eval]));
        
        P_j1_mag=P2_j1-P1_j1;
        P_j1_hub=P1_j1;
        
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
        C_10=double(subs(C10_temp,q, q_eval));
        
        C20_temp=X_1(2);
        C_20=double(subs(C20_temp,q, q_eval));
        
        D20_temp=X_1(3);
        D_20=double(subs(D20_temp,q, q_eval));
        
        C30_temp=X_1(4);
        C_30=double(subs(C30_temp,q, q_eval));
        
        D30_temp=X_1(5);
        D_30=double(subs(D30_temp,q, q_eval));
        
        C40_temp=X_1(6);
        C_40=double(subs(C40_temp,q, q_eval));
        
        D40_temp=X_1(7);
        D_40=double(subs(D40_temp,q, q_eval));
        
        
%         B_calculated=-q_eval*(C_20*R_2^(q_eval-1)+D_20*R_2^(-q_eval-1));
        B_calculated=-q_eval*(C_30*R_wave^(q_eval-1)+D_30*R_wave^(-q_eval-1));
        
        
        % CHANGES HERE
        K_B_ms=B_ms_given/B_calculated;
        
        
        %% Final estimation of the rotor losses:
        
        % Multiple slots:
        P_mag_ms=alpha_p*P_j1_mag*K_B_ms^2;
        P_hub_ms=alpha_p*P_j1_hub*K_B_ms^2;
        
        P_ms_matrix_n(count)=P_mag_ms;
        P_ms_tot_n=P_ms_tot+P_mag_ms;
        
        
        
        
    end




%% results


P_ms_tot
P_ms_matrix
vec_space2
vec_time2

P_ms_tot_n
P_ms_matrix_n
vec_space2_n
vec_time2_n

figure
hArray = bar(vec_space2,transpose(P_ms_matrix));
set(hArray(1),'FaceColor','w');

grid


ylabel('Rotor Loss (W)')
xlabel('Space order')
title('Temporal order 24')

legend('CP Multiple Slots')


