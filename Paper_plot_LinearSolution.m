%%%%%%%%%%%%%%%%%%%%%%%
%%% ENSO Complexity %%%
%%% Linear solution %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Solving the linear solution 
% Parameters
La = 8/3; % equatorial belt length
LaDim = 40000; % global
LoDim = 17500; % Pacific
Na = 64*2; % grid points for atmosphere model
No = round(LoDim/LaDim * Na); % grid points for ocean model
Lo = LoDim/LaDim*La;
dx = La/Na;
chi1 = 0.65/2; % projection from ocn to atm
chi2 = 0.65*2; % projection from atm to ocn
Qbar = 0.9; % mean vertical moisture gradient
c1 = 0.18; % take reduced gravity 0.03, mean thermocline depth 50, typical ocean velocity 0.5; Fr = U/sqrt(gh); 
gamma = 6.529; % wind stress coefficient
zeta = 8.7; %latent heating exchange coefficient
rW = 0.5; % Reflection at western boundary
rE = 1.0; % Reflection at eastern boundary

% Eq parameters
qc = 7; % latent heating multiplier coefficient
qe = 0.09296; % latent heating exponential coefficient
Tbar = 25/1.5; %1.5 is dimension, mean SST
tauq = 15; % latent heating adjustment rate
alpha_q = qc*qe*exp(qe*Tbar)/tauq; % latent heating factor
%alpha_q = 0.29;

% Thermocline coefficient eta depends on x
x = 0:dx:Lo-dx;
eta = 1.3 + 1.1 *tanh(7.5*(x-Lo/2))*0.95;  

% Zonal advective coefficient eta2 depends on x

%eta2 = max(0,4*exp(-(x-Lo/(7/3)).^2/0.1))*0.9*0.9;
eta2 = max(0,4*exp(-(x-Lo/(7/3)).^2/0.1))*0.9*1.3;
% Artificial damping
da = 1e-8;

% Redifined parameters
a = chi2 * c1 * gamma; 
b = c1 * zeta * alpha_q; 
m1 = - chi1 * alpha_q / (2 - 2 * Qbar);
m2 = - chi1 * alpha_q / (3 - 3 * Qbar);


% Linear matrix
M_11 = eye(No);
M_22 = -eye(No);
M_K = (- 1 + da * dx) * eye(Na);
M_R = (1 + 3 * da * dx) * eye(Na);


for i = 1:No-1
    M_11(i+1,i)= -1;
    M_22(i,i+1) = 1;
end
for i = 1:Na-1
    M_K(i,i+1) = 1;
    M_R(i,i+1) = -1;
end
M_11 = -c1 / dx * M_11;
M_22 = c1 / 3 / dx * M_22;

M_12 = zeros(No,No);
M_12(1,1) = c1/dx * rW;

M_21 = zeros(No,No);
M_21(end,end) = c1/3/dx * rE;

M_K(end,1) = 1;
M_R(end,1) = -1;

c1_tilde = (Na-1)/Na * dx * m1;
c2_tilde = 3*(Na-1)/Na * dx * m2;


B_K = -1/(Na-1) * ones(Na,No);


for i = 1:No
    B_K(i,i) = 1;
end



B_R = B_K;
M_13 = a/2 * (c1_tilde * inv(M_K) * B_K - c2_tilde * inv(M_R) * B_R);
M_23 = - a/3 * (c1_tilde * inv(M_K) * B_K - c2_tilde * inv(M_R) * B_R);

M_13 = M_13(1:No,:);
M_23 = M_23(1:No,:);

temp22 = (c1_tilde * inv(M_K) * B_K - c2_tilde * inv(M_R) * B_R);
M_33 = -b * eye(No);
M_31 = zeros(No,No);
M_32 = M_31;
for i = 1:No
    M_31(i,i) = c1*(eta(i)+eta2(i));
    M_32(i,i) = c1*(eta(i)-eta2(i));
end
% M_32 = M_31;

M = [M_11,M_12,M_13;
    M_21,M_22,M_23;
    M_31,M_32,M_33];

% For reconstruction the atmosphere variables
Ka_temp = c1_tilde * inv(M_K) * B_K;
Ra_temp = c2_tilde * inv(M_R) * B_R;

% Figure for the eigen values and eigen vectors
figure
[V, egvalue] = eig(M);
egvalue = diag(egvalue);
[freq,num] = sort(imag(egvalue),'ascend');
growth = real(egvalue(num));
[vlu, num2] = max(growth);
subplot(2,2,1)
hold on
plot(freq)
plot(num2,freq(num2),'o')
set(gca,'fontsize',12)
box on
title('Ranked frequency')
xlabel('Num')

subplot(2,2,2)
hold on
plot(growth)
plot(num2,growth(num2),'o')
xlabel('Num')
set(gca,'fontsize',12)
box on
title('Ranked growth rate')

[growth,num] = sort(real(egvalue),'descend');
freq = imag(egvalue(num));
subplot(2,2,3)
hold on
plot(freq(1:30)/(2*pi)*365/33,growth(1:30)*(365/33),'bo','linewidth',2)
plot(freq(1:2)/(2*pi)*365/33,growth(1:2)*(365/33),'ro','linewidth',2)
ylabel('Growth: yr^{-1}')
xlabel('Freq: yr^{-1}')
set(gca,'fontsize',12)
box on
title('First leading a few modes')

subplot(2,2,4)
hold on
xx1 = linspace(0,Lo,No);
ddx = xx1(2)-xx1(1);
xx2 = linspace(Lo+ddx,Lo*2,No);
xx3 = linspace(2*Lo+ddx,Lo*3,No);
plot(xx1,real(V(1:No,num(1))))
plot(xx1,imag(V(1:No,num(1))),'r')
plot(xx2,real(V(1+No:2*No,num(1))))
plot(xx2,imag(V(1+No:2*No,num(1))),'r')
plot(xx3,real(V(1+2*No:3*No,num(1))))
plot(xx3,imag(V(1+2*No:3*No,num(1))),'r')
legend('Real','Imag')
xlabel('x ([0,Lo], [0,Lo], [0,Lo])')
set(gca,'fontsize',12)
box on
title('Eigvector for enso mode 1')


figure
hold on
plot(x,eta,'b')
plot(x,eta2,'r')
% plot(x,eta3,'g')
legend('Thermocline', 'advection')

disp('Period of the leading mode:')
disp(2*pi/freq(1)*33/365)

%%%%%%%%%%%%%%%%%%%%% Meridional Basis
c = 0.05;
yy = -8:0.1:8;

% Atmopshere
phi_0 = 1/pi^(1/4) * exp(-yy.^2/2);
phi_2 = 1/pi^(1/4)/sqrt(8) * (4 * yy.^2 - 2) .* exp(-yy.^2/2);
% Ocean
psi_0 = 1/pi^(1/4) * exp(-yy.^2/2/c);
psi_2 = 1/pi^(1/4)/sqrt(8) * (4 * yy.^2/c - 2) .* exp(-yy.^2/2/c);

% Value at the equator
phi_0_eq = max(phi_0);
psi_0_eq = max(psi_0);

phi_2_eq = min(phi_2);
psi_2_eq = min(psi_2);

% Showing the meridional basis profiles
figure
subplot(2,2,1)
plot(yy*1.5,phi_0,'b','linewidth',2)
set(gca,'fontsize',12)
xlim([-7,7])
box on
title('Atm \phi_0')
subplot(2,2,2)
plot(yy*1.5,phi_2,'b','linewidth',2)
set(gca,'fontsize',12)
xlim([-7,7])
box on
title('Atm \phi_2')
subplot(2,2,3)
plot(yy*1.5,psi_0,'b','linewidth',2)
set(gca,'fontsize',12)
box on
xlim([-7,7])
title('Ocn \psi_0')
subplot(2,2,4)
plot(yy*1.5,psi_2,'b','linewidth',2)
set(gca,'fontsize',12)
box on
xlim([-7,7])
title('Ocn \psi_2')


%%%%%%%%%%%%%%%%% Spatial reconstruction of the linear solutions
% Note that the decaying rate is set to be zero in these figures for the
% convenience of presenting the results
 
% Units
Dim_T = 1.5;
Dim_U = 0.5;
Dim_H = 20.8;
Dim_u = 50;
Dim_x = 40000/Na;

% Four Enso bases
vec1_Ko = V(1:No,num(1));
vec1_Ro = V(No+1:2*No,num(1));
vec1_T = V(2*No+1:3*No,num(1));

vec2_Ko = V(1:No,num(2));
vec2_Ro = V(No+1:2*No,num(2));
vec2_T = V(2*No+1:3*No,num(2));

 
% Time step for figure not for numerical intergration
dt = (1/freq(1))/100*1.5;
T = 68;

dn = 1;
LL = round(T/dt/dn);
Hov_Ko = zeros(LL,No);
Hov_Ro = zeros(LL,No);
Hov_Ka = zeros(LL,Na);
Hov_Ra = zeros(LL,Na);
Hov_T = zeros(LL,No);
i = 1;
for t = dt:dn*dt:dn*dt*LL
    Hov_T(i,:) = 10/2 * Dim_T * psi_0_eq * (exp(1i * freq(1) * t) * vec1_T...
        +  exp( 1i*freq(2) * t) * vec2_T);
    Hov_Ko(i,:) = 10/2 * (exp( 1i * freq(1) * t) * vec1_Ko  +  exp( 1i * freq(2) * t)...
        * vec2_Ko);
    Hov_Ro(i,:) = 10/2 * (exp( 1i * freq(1) * t) * vec1_Ro  +  exp( 1i * freq(2) * t)...
        * vec2_Ro);
    Hov_Ka(i,:) = 10/2 * (Ka_temp * exp(1i * freq(1) * t) * vec1_T + Ka_temp *...
        exp(1i * freq(2) * t) * vec2_T);
    Hov_Ra(i,:) = 10/2 * (Ra_temp * exp(1i * freq(1) * t) * vec1_T +...
        Ra_temp * exp(1i * freq(2) * t) * vec2_T);
    i = i + 1;
end
Hov_U = Dim_U* ((Hov_Ko - Hov_Ro) * psi_0_eq + Hov_Ro/sqrt(2) * psi_2_eq);
Hov_H = Dim_H* ((Hov_Ko + Hov_Ro) * psi_0_eq + Hov_Ro/sqrt(2) * psi_2_eq);
Hov_u = Dim_u* ((Hov_Ka - Hov_Ra) * phi_0_eq + Hov_Ra/sqrt(2) * phi_2_eq);
Hov_u = real(Hov_u(:,1:No));




figure 
colormap('jet');
subplot(1,4,1)
levels=-20:.2:20;
contourf( (0:No-1)*625/1000/2, (dt:dn*dt:dn*dt*LL)*34/365, Hov_u,'linestyle','none')
colorbar('southoutside')
title('u','fontsize',16);
set(gca,'fontsize',12)
ylabel('Year')
set(gca,'XTick',[1,21,41]*625/1000/2);
set(gca,'XTicklabel',{'120','180','240'});
subplot(1,4,2)
levels=-2:.1:2;
contourf( (0:No-1)*625/1000/2, (dt:dn*dt:dn*dt*LL)*34/365, Hov_U,'linestyle','none')
colorbar('southoutside')
title('U','fontsize',16);
set(gca,'fontsize',12)
ylabel('Year')
set(gca,'XTick',[1,21,41]*625/1000/2);
set(gca,'XTicklabel',{'120','180','240'});
subplot(1,4,3)
levels=-60:1:60;
contourf( (0:No-1)*625/1000/2, (dt:dn*dt:dn*dt*LL)*34/365, Hov_H,'linestyle','none')
colorbar('southoutside')
title('H','fontsize',16);
set(gca,'fontsize',12)
ylabel('Year')
set(gca,'XTick',[1,21,41]*625/1000/2);
set(gca,'XTicklabel',{'120','180','240'});
subplot(1,4,4)
levels=-10:1:10;
contourf( (0:No-1)*625/1000/2, (dt:dn*dt:dn*dt*LL)*34/365, Hov_T,'linestyle','none')
colorbar('southoutside')
title('T','fontsize',16);
set(gca,'fontsize',12)
ylabel('Year')
set(gca,'XTick',[1,21,41]*625/1000/2);
set(gca,'XTicklabel',{'120','180','240'});
%sgtitle('(a) EP El Nino Dominate Mode')
sgtitle('(b) CP El Nino Dominate Mode')