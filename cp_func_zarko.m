function [ output_args ] = cp_func_zarko( Rm,Rs,Rr,bo,Qs,R )
%Calculates the CP according to Zarko.
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
alphao=asin(bo/Rs); %slot opening in rad
theta1=thetas/2-alphao/2;
theta2=theta1+alphao;
%----------------------------------------------------------------------
%CONFORMAL MAPPING OF THE SLOT STRUCTURE
clear j
%Coefficients of the conformal mapping function
boprime=theta2-theta1;
gprime=log(Rp2/Rp1);
bCM=(boprime/2/gprime+sqrt((boprime/2/gprime)^2+1))^2;
aCM=1/bCM;
%CALCULATE COMPLEX RELATIVE AIR GAP PERMEANCE
clear z_real z_imag
z_imag=j*[0:thetas/N:thetas];
z_real(1:length(z_imag))=log(R);
z=z_real+z_imag; %Coordinates of the points in the Z plane
s=exp(z);
%Define conformal mapping function
FCM=inline('j*gprime/pi*(log((1+s)/(1-s))-log((b+s)/(b-s))-2*(b-1)/sqrt(b)*atan(s/sqrt(b)))+log(Rp2)+j*theta2',...
    's','b','gprime','Rp2','theta2');
%START OPTIMIZATION
w0=[0.01*aCM]; %initial values
for i=1:length(z)
    res=@(w0)(FCM(sqrt((w0-bCM)/(w0-aCM)),bCM,gprime,Rp2,theta2)-z(i));
    [w(i),resnorm,residual,exitflag]=lsqnonlin(res,w0,[],[]); % Call optimizer
    k(i)=exp(j*gprime/pi*log(w(i))+log(Rp2)+j*thetas/2);
    lambda(i)=k(i)*(w(i)-1)/(w(i)-aCM)^(1/2)/(w(i)-bCM)^(1/2)/s(i);
    %fprintf(1,'Iteration: %d, dt/dz: %f\n',i,lambda(i));
    w0=w(i);
end
lambda=[conj(lambda(N/2+1:-1:1)) lambda(2:N/2)];
angle=(1:length(real(lambda)))*360/length(real(lambda))-180;
%figure;plot(angle,real(lambda));grid
%figure;plot(angle,imag(lambda));grid

clear k
for k=2:length(lambda)
    if imag(lambda(k))<0 && imag(lambda(k+1))>0
        flag=k;
    end
end

for k=1:flag
    lambda_zarko(2*flag-k+1)=lambda(flag-k+1);
    lambda_zarko(k+1)=real(lambda(flag-k+1))-j*imag(lambda(flag-k+1));
    
end
lambda_zarko(1)=conj(lambda_zarko(length(lambda_zarko)));
angle_zarko=(0:(2*flag-1))*360/2/flag;


output_args=[angle_zarko; lambda_zarko];


end

