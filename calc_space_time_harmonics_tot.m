
% Padova machine:
% Calculation of the space and time harmonics.

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
N_points=1000; % number of points for the waveform

r=R_m+.05*10^-3; % [m] radius where we want to calculate the waveform (used in 2DFEA)
%r=28.502*10^-3; % [m] radius where we want to calculate the waveform (used in 2DFEA)
R_wave=r;




N_pos=555;



%% Calculation of the coefficientsof magnetisation:

M_vec=zeros(1,N_harm);

for count=1:2:N_harm
    n=count;
    A1(count)=sin((n*p+1)*alpha_p*pi/2/p)/((n*p+1)*alpha_p*pi/2/p);
    A2(count)=sin((n*p-1)*alpha_p*pi/2/p)/((n*p-1)*alpha_p*pi/2/p);
    
    Mr_vec(count)=(B_r/mu_0)*alpha_p*(A1(count)+A2(count));
    Mt_vec(count)=(B_r/mu_0)*alpha_p*(A1(count)-A2(count));
    
    M_vec(count)=Mr_vec(count)+n*p*Mt_vec(count);
    
    A3(count)=(n*p-1/n/p)*Mr_vec(count)/M_vec(count)+1/n/p;
end




%% Calculate waveform for a given r>Rm:
% In the air-gap.

if R_m<r
    % eq. (34) and (35) in Zhu's paper 1993:
    
    theta=(0:N_points)*pi/p/(N_points)-pi/2/p;
    %theta=(0:N_points)*2*pi/(N_points+1); % theta in mech radians
    
    % Calculation of the coeficients:
    
    for count=1:2:N_harm
        n=count;
        KB=mu_0*M_vec(count)/mu_r*n*p/((n*p)^2-1)*((A3(count)-1)+2*(R_r/R_m)^(n*p+1)-(A3(count)+1)*(R_r/R_m)^(2*n*p))/((mu_r+1)/mu_r*(1-(R_r/R_s)^(2*n*p))-(mu_r-1)/mu_r*((R_m/R_s)^(2*n*p)-(R_r/R_m)^(2*n*p))); % checked
        f_Br=(r/R_s)^(n*p-1)*(R_m/R_s)^(n*p+1)+(R_m/r)^(n*p+1);
        f_Bt=-(r/R_s)^(n*p-1)*(R_m/R_s)^(n*p+1)+(R_m/r)^(n*p+1);
        
        coefs_radial(count)=KB*f_Br;
        coefs_tan(count)=KB*f_Bt;
        
%         B_rI_temp(count,:)=KB*f_Br*cos(n*p*theta);
%         B_tI_temp(count,:)=-KB*f_Br*sin(n*p*theta);
    end
    
%     figure
%     plot(coefs_radial,'o')
    
    
    
    B_rI=zeros(1,length(theta));
    B_tI=zeros(1,length(theta));
    for count=1:N_harm
        n=count;
        B_rI=B_rI+coefs_radial(count)*cos(n*p*theta);
        B_tI=B_tI+coefs_tan(count)*sin(n*p*theta);
        
    end
    B_r=B_rI;
    B_t=B_tI;
    
end

%% Calculate waveform for a given r<=Rm
% In the magnets. FINISH THIS LATER

if R_m>=r
    
    % eq. (36) and (37) in Zhu's paper 1993:
    
    theta=(0:N_points)*pi/p/(N_points)+pi/2/p/2;
    %theta=(0:N_points)*2*pi/(N_points+1); % theta in mech radians
    
    % Calculation of the coeficients:
    
    for count=1:2:N_harm
        n=count;
        C2=((A3(count)-1/mu_r)*(R_m/R_s)^(2*n*p)+(1+1/mu_r)*(R_r/R_m)^(n*p+1)*(R_m/R_s)^(2*n*p)-(A3(count)+1/mu_r)-(1-1/mu_r)*(R_r/R_m)^(n*p+1))/((mu_r+1)/mu_r*(1-(R_r/R_s)^(2*n*p))-(mu_r-1)/mu_r*((R_m/R_s)^(2*n*p)-(R_r/R_m)^(2*n*p)));
        
        C1=mu_0*M_vec(count)*n*p/((n*p)^2-1)*C2*((r/R_m)^(n*p-1)+(R_r/R_m)^(n*p-1)*(R_r/r)^(n*p+1));
        
        C1s=mu_0*M_vec(count)*n*p/((n*p)^2-1)*C2*((r/R_m)^(n*p-1)-(R_r/R_m)^(n*p-1)*(R_r/r)^(n*p+1));
        
        B1=mu_0*M_vec(count)*n*p/((n*p)^2-1)*(R_r/r)^(n*p+1);
        
        D1=mu_0*M_vec(count)*n*p/((n*p)^2-1)*A3(count);
        D1s=mu_0*M_vec(count)/((n*p)^2-1)*A3(count);
        
        coefs_radial(count)=C1+B1+D1;
        coefs_tan(count)=-C1s+B1-D1s;
        
%         B_rI_temp(count,:)=KB*f_Br*cos(n*p*theta);
%         B_tI_temp(count,:)=-KB*f_Br*sin(n*p*theta);
    end
    
%     figure
%     plot(coefs_radial,'o')
    
    
    
    B_rII=zeros(1,length(theta));
    B_tII=zeros(1,length(theta));
    for count=1:N_harm
        n=count;
        B_rII=B_rII+coefs_radial(count)*cos(n*p*theta);
        B_tII=B_tII+coefs_tan(count)*sin(n*p*theta);
        
    end
    
    B_r=B_rII;
    B_t=B_tII;
    
end


theta_p=(0:N_points)*2*pi/p/(N_points); % electrical angle




Br_slotless=B_r;
Bt_slotless=B_t;
theta=theta_p; % electrical angle


%% Complex permeance:







temp2=cp_func_renedo( R_m,R_s,R_r,bo,Qs,R_wave,dslot );
cp_func=temp2(2,:);
theta_points=temp2(1,:); % mechanical angle


% Calculate the coeffs of the fourier series multiple slots:

f1_t=imag(cp_func);
f1_r=real(cp_func);

theta_integration=theta_points*Qs; % to make it up to 2*pi

a_0=1/2/pi*trapz(theta_integration,f1_r);
a_1=1/pi*trapz(theta_integration,f1_r.*cos(1*theta_integration));
a_2=1/pi*trapz(theta_integration,f1_r.*cos(2*theta_integration));
a_3=1/pi*trapz(theta_integration,f1_r.*cos(3*theta_integration));
a_4=1/pi*trapz(theta_integration,f1_r.*cos(4*theta_integration));
a_5=1/pi*trapz(theta_integration,f1_r.*cos(5*theta_integration));
a_6=1/pi*trapz(theta_integration,f1_r.*cos(6*theta_integration));
a_7=1/pi*trapz(theta_integration,f1_r.*cos(7*theta_integration));

b_1=1/pi*trapz(theta_integration,f1_t.*sin(1*theta_integration));
b_2=1/pi*trapz(theta_integration,f1_t.*sin(2*theta_integration));
b_3=1/pi*trapz(theta_integration,f1_t.*sin(3*theta_integration));
b_4=1/pi*trapz(theta_integration,f1_t.*sin(4*theta_integration));
b_5=1/pi*trapz(theta_integration,f1_t.*sin(5*theta_integration));
b_6=1/pi*trapz(theta_integration,f1_t.*sin(6*theta_integration));
b_7=1/pi*trapz(theta_integration,f1_t.*sin(7*theta_integration));

% make sure the phase is right.
alpha=0;
theta_temp=Qs*(theta/p+alpha);
f_tempr=a_0+a_1*cos(theta_temp)+a_2*cos(2*theta_temp)+a_3*cos(3*theta_temp)+a_4*cos(4*theta_temp)+a_5*cos(5*theta_temp)+a_6*cos(6*theta_temp)+a_7*cos(7*theta_temp);
f_tempt=b_1*sin(theta_temp)+b_2*sin(2*theta_temp)+b_3*sin(3*theta_temp)+b_4*sin(4*theta_temp)+b_5*sin(5*theta_temp)+b_6*sin(6*theta_temp)+b_7*sin(7*theta_temp);

% Apply the complex permeance:

B_slotless_complex=Br_slotless+i*Bt_slotless;
cp_func_adapted=conj(f_tempr+i*f_tempt);





Bcp_complex=B_slotless_complex.*cp_func_adapted;
Bcp_r=real(Bcp_complex);
Bcp_t=imag(Bcp_complex);




%% Several rotor positions:

clear theta theta_temp
% make sure the phase is right.

%  Multiple slots:




alpha_vec=(0:N_pos)/(N_pos)*2*pi/Qs; % mechanical

theta=(0:N_points)*pi/(N_points); % electrical
alpha=0;


count=1;




x_ms=zeros(length(theta),length(alpha_vec));
y_ms=zeros(length(theta),length(alpha_vec));


for alpha=alpha_vec;
    theta_temp=Qs*(theta/p+alpha); % to make it up to 2*pi
    f_tempr=a_0+a_1*cos(theta_temp)+a_2*cos(2*theta_temp)+a_3*cos(3*theta_temp)+a_4*cos(4*theta_temp)+a_5*cos(5*theta_temp)+a_6*cos(6*theta_temp)+a_7*cos(7*theta_temp);
    f_tempt=b_1*sin(theta_temp)+b_2*sin(2*theta_temp)+b_3*sin(3*theta_temp)+b_4*sin(4*theta_temp)+b_5*sin(5*theta_temp)+b_6*sin(6*theta_temp)+b_7*sin(7*theta_temp);
    
    
    
    
    % Apply the complex permeance:
    
    B_slotless_complex=Br_slotless+i*Bt_slotless;
    cp_func_adapted=conj(f_tempr+i*f_tempt);
    

    
    
    
    Bcp_complex=B_slotless_complex.*cp_func_adapted;
    Bcp_r=real(Bcp_complex);
    Bcp_t=imag(Bcp_complex);
    
    x_ms(:,count)=Bcp_r;
    y_ms(:,count)=Bcp_t;
    
    
    
    count=count+1;
    
end






time=transpose(theta/pi*180);
count=1;
for alpha=alpha_vec;
    
    count;
    
    count=count+1;
end



%% Calculate FFT



% Multiple slots:

[m1,n1]=size(x_ms);
%Matrix manipulation
%Data for N-pole Xn
%
x1=(x_ms(1,:)-x_ms(m1,:))/2;
x1n=x_ms(2:m1-1,:);
Xn=[x1;x1n];
%Generating points for S-pole, data for S-pole is the negative of data for N-pole
Xs=-Xn;
%
Xns=[Xn;Xs];
time2 = linspace(0,time(end)*2,size(Xns,1));

[m,n]=size(Xns);
Xnsfft2=fft2(Xns);
Xnsabs=abs(Xnsfft2)/(m*n);
%Matrix manipulation
X=Xnsabs;
if rem(m,2)==0,
    X1m0=X(2:(m/2),1);
    X2m0=X(m/2+2:m,1);
    X2m0=flipud(X2m0);
    Xm0=X1m0+X2m0;
    %
    X13=X(2:m/2,2:n);
    X24=X(m/2+2:m,2:n);
    X24=flipud(fliplr(X24));
    X1234=X13+X24;
    Y=zeros(m/2,n);
    Y(1,:)=X(1,:);
    Y(2:m/2,1)=Xm0;
    Y(2:m/2,2:n)=X1234;
else
    X1m0=X(2:m/2+0.5,1);
    X2m0=X(m/2+1:m,1);
    X2m0=flipud(X2m0);
    Xm0=X1m0+X2m0;
    %
    X13=X(2:m/2+0.5,2:n);
    X24=X(m/2+1:m,2:n);
    X24=flipud(fliplr(X24));
    X1234=X13+X24;
    %
    Y=zeros(m/2+0.5,n);
    Y(1,:)=X(1,:);
    Y(2:m/2+0.5,1)=Xm0;
    Y(2:m/2+0.5,2:n)=X1234;
end
% if Y is partitioned from the middle column, then the first part represent the forward rotating waves and second part the backward rotating waves.
Xfft2=Y;
%Mat_to_write = Y(2:2:20,1:6);
Mat_to_write = Y;
%xlswrite('multiple_slots_mat',Mat_to_write)

Mat_multiple_slots=Y(1:70,1:70);

for k1=1:length(Y(:,1))
    %Y_n(k1,:)=Y(k1,:);
    for k2=1:length(Y(1,:))
        Y_n(k1,k2)=Y(k1,length(Y(1,:))-k2+1);
    end
end
Mat_multiple_slots_n=Y_n(1:70,1:70);


%% Complex permeance: NO CURVATURE





temp2=cp_func_renedo_no_curv( R_m,R_s,R_r,bo,Qs,R_wave,dslot );
cp_func=temp2(2,:);
theta_points=temp2(1,:); % mechanical angle


% Calculate the coeffs of the fourier series multiple slots:

f1_t=imag(cp_func);
f1_r=real(cp_func);

theta_integration=theta_points*Qs; % to make it up to 2*pi

a_0=1/2/pi*trapz(theta_integration,f1_r);
a_1=1/pi*trapz(theta_integration,f1_r.*cos(1*theta_integration));
a_2=1/pi*trapz(theta_integration,f1_r.*cos(2*theta_integration));
a_3=1/pi*trapz(theta_integration,f1_r.*cos(3*theta_integration));
a_4=1/pi*trapz(theta_integration,f1_r.*cos(4*theta_integration));
a_5=1/pi*trapz(theta_integration,f1_r.*cos(5*theta_integration));
a_6=1/pi*trapz(theta_integration,f1_r.*cos(6*theta_integration));
a_7=1/pi*trapz(theta_integration,f1_r.*cos(7*theta_integration));

b_1=1/pi*trapz(theta_integration,f1_t.*sin(1*theta_integration));
b_2=1/pi*trapz(theta_integration,f1_t.*sin(2*theta_integration));
b_3=1/pi*trapz(theta_integration,f1_t.*sin(3*theta_integration));
b_4=1/pi*trapz(theta_integration,f1_t.*sin(4*theta_integration));
b_5=1/pi*trapz(theta_integration,f1_t.*sin(5*theta_integration));
b_6=1/pi*trapz(theta_integration,f1_t.*sin(6*theta_integration));
b_7=1/pi*trapz(theta_integration,f1_t.*sin(7*theta_integration));

% make sure the phase is right.
alpha=0;
theta_temp=Qs*(theta/p+alpha);
f_tempr=a_0+a_1*cos(theta_temp)+a_2*cos(2*theta_temp)+a_3*cos(3*theta_temp)+a_4*cos(4*theta_temp)+a_5*cos(5*theta_temp)+a_6*cos(6*theta_temp)+a_7*cos(7*theta_temp);
f_tempt=b_1*sin(theta_temp)+b_2*sin(2*theta_temp)+b_3*sin(3*theta_temp)+b_4*sin(4*theta_temp)+b_5*sin(5*theta_temp)+b_6*sin(6*theta_temp)+b_7*sin(7*theta_temp);

% Apply the complex permeance:

B_slotless_complex=Br_slotless+i*Bt_slotless;
cp_func_adapted=conj(f_tempr+i*f_tempt);





Bcp_complex=B_slotless_complex.*cp_func_adapted;
Bcp_r=real(Bcp_complex);
Bcp_t=imag(Bcp_complex);




%% Several rotor positions:

clear theta theta_temp
% make sure the phase is right.

%  Multiple slots:




alpha_vec=(0:N_pos)/(N_pos)*2*pi/Qs; % mechanical

theta=(0:N_points)*pi/(N_points); % electrical
alpha=0;


count=1;




x_nc=zeros(length(theta),length(alpha_vec));
y_nc=zeros(length(theta),length(alpha_vec));


for alpha=alpha_vec;
    theta_temp=Qs*(theta/p+alpha); % to make it up to 2*pi
    f_tempr=a_0+a_1*cos(theta_temp)+a_2*cos(2*theta_temp)+a_3*cos(3*theta_temp)+a_4*cos(4*theta_temp)+a_5*cos(5*theta_temp)+a_6*cos(6*theta_temp)+a_7*cos(7*theta_temp);
    f_tempt=b_1*sin(theta_temp)+b_2*sin(2*theta_temp)+b_3*sin(3*theta_temp)+b_4*sin(4*theta_temp)+b_5*sin(5*theta_temp)+b_6*sin(6*theta_temp)+b_7*sin(7*theta_temp);
    
    
    
    
    % Apply the complex permeance:
    
    B_slotless_complex=Br_slotless+i*Bt_slotless;
    cp_func_adapted=conj(f_tempr+i*f_tempt);
    
    
    
    
    Bcp_complex=B_slotless_complex.*cp_func_adapted;
    Bcp_r=real(Bcp_complex);
    Bcp_t=imag(Bcp_complex);
    
    x_nc(:,count)=Bcp_r;
    y_nc(:,count)=Bcp_t;
    
    
    
    count=count+1;
    
end






time=transpose(theta/pi*180);
count=1;
for alpha=alpha_vec;
    
    count;
    
    count=count+1;
end



%% Calculate FFT



% Multiple slots:

[m1,n1]=size(x_nc);
%Matrix manipulation
%Data for N-pole Xn
%
x1=(x_nc(1,:)-x_nc(m1,:))/2;
x1n=x_nc(2:m1-1,:);
Xn=[x1;x1n];
%Generating points for S-pole, data for S-pole is the negative of data for N-pole
Xs=-Xn;
%
Xns=[Xn;Xs];
time2 = linspace(0,time(end)*2,size(Xns,1));

[m,n]=size(Xns);
Xnsfft2=fft2(Xns);
Xnsabs=abs(Xnsfft2)/(m*n);
%Matrix manipulation
X=Xnsabs;
if rem(m,2)==0,
    X1m0=X(2:(m/2),1);
    X2m0=X(m/2+2:m,1);
    X2m0=flipud(X2m0);
    Xm0=X1m0+X2m0;
    %
    X13=X(2:m/2,2:n);
    X24=X(m/2+2:m,2:n);
    X24=flipud(fliplr(X24));
    X1234=X13+X24;
    Y=zeros(m/2,n);
    Y(1,:)=X(1,:);
    Y(2:m/2,1)=Xm0;
    Y(2:m/2,2:n)=X1234;
else
    X1m0=X(2:m/2+0.5,1);
    X2m0=X(m/2+1:m,1);
    X2m0=flipud(X2m0);
    Xm0=X1m0+X2m0;
    %
    X13=X(2:m/2+0.5,2:n);
    X24=X(m/2+1:m,2:n);
    X24=flipud(fliplr(X24));
    X1234=X13+X24;
    %
    Y=zeros(m/2+0.5,n);
    Y(1,:)=X(1,:);
    Y(2:m/2+0.5,1)=Xm0;
    Y(2:m/2+0.5,2:n)=X1234;
end
% if Y is partitioned from the middle column, then the first part represent the forward rotating waves and second part the backward rotating waves.
Xfft2=Y;
%Mat_to_write = Y(2:2:20,1:6);
Mat_to_write = Y;
%xlswrite('multiple_slots_mat',Mat_to_write)

Mat_no_curv=Y(1:70,1:70);

for k1=1:length(Y(:,1))
    %Y_n(k1,:)=Y(k1,:);
    for k2=1:length(Y(1,:))
        Y_n(k1,k2)=Y(k1,length(Y(1,:))-k2+1);
    end
end
Mat_no_curv_n=Y_n(1:70,1:70);


%% Save the results in matrices:
% careful with the harmonics!

vec_time=[0 Qs/p 2*Qs/p 3*Qs/p];
vec_space=1:2:57;

length2=2*length(vec_space);

Mat_multiple_slot_export=zeros(length2/2,3);
Mat_multiple_slot_export(:,1)=Mat_multiple_slots(2:2:length2,1);
Mat_multiple_slot_export(:,2)=Mat_multiple_slots(2:2:length2,2);
Mat_multiple_slot_export(:,3)=Mat_multiple_slots(2:2:length2,3);
Mat_multiple_slot_export(:,4)=Mat_multiple_slots(2:2:length2,4);

Mat_n=zeros(length2/2,3);
Mat_no_curv_export(:,1)=Mat_no_curv(2:2:length2,1);
Mat_no_curv_export(:,2)=Mat_no_curv(2:2:length2,2);
Mat_no_curv_export(:,3)=Mat_no_curv(2:2:length2,3);
Mat_no_curv_export(:,4)=Mat_no_curv(2:2:length2,4);



%% Save results


save('data_to_plot','vec_time','vec_space','Mat_multiple_slot_export','Mat_no_curv_export')

