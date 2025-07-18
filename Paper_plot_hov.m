EigenSolver_rescale

Ua_store = []; Ka_bar_store = []; Ra_bar_store = []; Ko_store = []; Ro_store = []; A_store = []; T_store = []; Q_store = []; A_phy_store = []; K_phy_store = []; ap_store = [];
R_phy_store = []; I_store = []; A_bar_store =[]; qq_store = [];
noise_record=[];

num = 200;

for yy = 1:num
    filename_ENSO = ['Simulations\main_test5new_mu_',num2str(yy),'.mat'];
    load(filename_ENSO)
    Ua_store = [Ua_store,Ua];
    A_store = [A_store,A];
    A_bar_store = [A_bar_store,A_bar_record];
    Ko_store = [Ko_store,Ko];
    Ro_store = [Ro_store,Ro];
    T_store = [T_store,T];
    Q_store = [Q_store,Q_phy];
    A_phy_store = [A_phy_store,A_phy];
    K_phy_store = [K_phy_store,K_phy];
    R_phy_store = [R_phy_store,R_phy];
    I_store = [I_store, I_record];
    noise_record = [noise_record,noise_temp_record];
    ap_store = [ap_store,ap];
    clear Ko Ro A T Ua Q_phy A_phy K_phy R_phy I_record A_bar_record noise_temp_record ap %Eq_record 
end

factor = 1;


%Ua = Ua_store;
Ko = Ko_store; Ro = Ro_store; T = T_store; Q_phy = Q_store;
K_phy = K_phy_store; R_phy = R_phy_store; A_phy = A_phy_store; I = I_store;
Ua = Ua_store; Ka_bar = Ka_bar_store; Ra_bar = Ra_bar_store; 
A_bar = A_bar_store;
ap_record = ap_store;

K = fftshift(fft(K_phy, [], 1), 1); R = fftshift(fft(R_phy, [], 1), 1); Q = fftshift(fft(Q_phy, [], 1), 1); A = fftshift(fft(A_phy, [], 1), 1);


dim_t = 34;
La = 8/3;
LaDim = 40000;
LoDim = 17500;
dn = 2; % a factor that can refine the grids
Na = 64 * dn; % atmosphere grid points
No = round(LoDim/LaDim*Na);
Lo = LoDim/LaDim*La;
dt = 8/24/dim_t/10; % time step, every dt*dim_t hours  
dx = La/Na; % distance between every two grid points
x = (0:Na-1) * dx; % grid points in the x axis
da = 10e-8;
dim_u = 5; % dimensional unit of velocity
dim_Q = 1.5; % dimensional unit of Q
dim_theta = 1.5; % dimensional unit of theta
dim_A = 1.5; % dimensional unit of A
Dim_T = 1.5;
Dim_H = 20.8;
Dim_U = 0.5;
gap = 30;


% Eq parameters
qc = 7; % latent heating multiplier coefficient
qe = 0.09296; % latent heating exponential coefficient
Tbar = 25/1.5; %1.5 is dimension, mean SST
tauq = 15; % latent heating adjustment rate
alpha_q = qc*qe*exp(qe*Tbar)/tauq; % latent heating factor

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

Hbar = 22;

dim_x = 15000; % one dimension unit of x axis (km) 



%dim_u = 50;
A_bar_Fourier = real(fftshift(fft(A_bar)));
u_bar = ( (Ua)*psi_0_eq) * dim_u; % unit 50m/s
A = A * dim_A * psi_0_eq; % unit 10K^-1
Hov_T_a = Dim_T * psi_0_eq * T;
Hov_U_a = Dim_U * ((Ko - Ro) * psi_0_eq + Ro/sqrt(2) * psi_2_eq);
Hov_H_a = Dim_H * ((Ko + Ro) * psi_0_eq + Ro/sqrt(2) * psi_2_eq);

T_E = mean(Hov_T_a(No/2+1:No,:));


u_phy = ( (K_phy-R_phy)*psi_0_eq + 1/sqrt(2)*R_phy*psi_2_eq ) * dim_u; % unit 50m/s
theta_phy = ((-K_phy-R_phy)*psi_0_eq - 1/sqrt(2)*R_phy*psi_2_eq) * dim_theta; % unit 15K 

Q_phy = Q_phy * dim_Q * psi_0_eq; % unit 15K
A_phy = A_phy * dim_A * psi_0_eq; % unit 10K^-1
A_bar = A_bar * dim_A * psi_0_eq; % unit 10K^-1


jj = size(K,2);

if mod(jj,2) == 0
    zero_mode2 = (jj)/2+1;
else
    zero_mode2 = (jj+1)/2;
end

if mod(Na,2) == 0
    zero_mode = Na/2+1;
else
    zero_mode = (Na+1)/2;
end




MJO1 = zeros(Na, jj); MJO2 = zeros(Na, jj);
MJO1_temp = zeros(Na, jj); MJO2_temp = zeros(Na, jj);

for j = 1:jj
    MJO1_Fourier = zeros(Na, 1);
    MJO2_Fourier = zeros(Na, 1);

    MJO1_Fourier(zero_mode+1) =  K(zero_mode+1,j) * eg_K_MJO(1)' + R(zero_mode+1,j) * eg_R_MJO(1)' + Q(zero_mode+1,j) * eg_Q_MJO(1)' + (A(zero_mode+1,j)+0*A_bar_Fourier(zero_mode+1)) * eg_A_MJO(1)';
    MJO1_Fourier(zero_mode+2) =  K(zero_mode+2,j) * eg_K_MJO(2)' + R(zero_mode+2,j) * eg_R_MJO(2)' + Q(zero_mode+2,j) * eg_Q_MJO(2)' + (A(zero_mode+2,j)+0*A_bar_Fourier(zero_mode+2)) * eg_A_MJO(2)';
    MJO1_Fourier(zero_mode+3) =  K(zero_mode+3,j) * eg_K_MJO(3)' + R(zero_mode+3,j) * eg_R_MJO(3)' + Q(zero_mode+3,j) * eg_Q_MJO(3)' + (A(zero_mode+3,j)+0*A_bar_Fourier(zero_mode+3)) * eg_A_MJO(3)';
    MJO1_temp(:,j) = ifft(ifftshift(MJO1_Fourier)); 

    MJO2_Fourier(zero_mode-1) =  K(zero_mode-1,j) * eg_K_MJO(1) + R(zero_mode-1,j) * eg_R_MJO(1) + Q(zero_mode-1,j) * eg_Q_MJO(1) + (A(zero_mode-1,j)+0*A_bar_Fourier(zero_mode+1)) * eg_A_MJO(1);
    MJO2_Fourier(zero_mode-2) =  K(zero_mode-2,j) * eg_K_MJO(2) + R(zero_mode-2,j) * eg_R_MJO(2) + Q(zero_mode-2,j) * eg_Q_MJO(2) + (A(zero_mode-2,j)+0*A_bar_Fourier(zero_mode+2)) * eg_A_MJO(2);
    MJO2_Fourier(zero_mode-3) =  K(zero_mode-3,j) * eg_K_MJO(3) + R(zero_mode-3,j) * eg_R_MJO(3) + Q(zero_mode-3,j) * eg_Q_MJO(3) + (A(zero_mode-3,j)+0*A_bar_Fourier(zero_mode+3)) * eg_A_MJO(3);
    MJO2_temp(:,j) = ifft(ifftshift(MJO2_Fourier));
end

bd1 = -90/(dim_t*dt*gap); bd2 = -30/(dim_t*dt*gap);
for i = 1:Na
    MJO1_temporal = fftshift(fft(MJO1_temp(i,:)));
    MJO1_temporal2 = MJO1_temporal * 0;
    MJO1_temporal2(zero_mode2+round(jj/bd2): zero_mode2+round(jj/bd1)) = MJO1_temporal(zero_mode2+round(jj/bd2): zero_mode2+round(jj/bd1));
    MJO1(i,:) = real(ifft(ifftshift(MJO1_temporal2)));
end

bd1 = 30/(dim_t*dt*gap); bd2 = 90/(dim_t*dt*gap);
for i = 1:Na
    MJO2_temporal = fftshift(fft(MJO2_temp(i,:)));
    MJO2_temporal2 = MJO2_temporal * 0;
    MJO2_temporal2(zero_mode2+round(jj/bd2): zero_mode2+round(jj/bd1)) = MJO2_temporal(zero_mode2+round(jj/bd2): zero_mode2+round(jj/bd1));
    MJO2(i,:) = real(ifft(ifftshift(MJO2_temporal2)));
end
MJO = MJO1+MJO2;

if floor(num/10) < 2
    listn = 0;
else
    listn = 0:floor(num/10)-1;
end


jj = size(K,2);



if mod(jj,2) == 0
    zero_mode2 = (jj)/2+1;
else
    zero_mode2 = (jj+1)/2;
end

if mod(Na,2) == 0
    zero_mode = Na/2+1;
else
    zero_mode = (Na+1)/2;
end



west_node = 1:round(No*1/2);
east_node = round(No*1/2)+1:No;
T3_node = 31:53;
T4_node = 14:31;
T34_node = 25:43;


T4_store_a = mean(T(T4_node,:)) * Dim_T * psi_0_eq;
Tw_store_a = mean(T(west_node,:)) * Dim_T * psi_0_eq;
T34_store_a = mean(T(T34_node,:)) * Dim_T * psi_0_eq;
Te_store_a = mean(T(east_node,:)) * Dim_T * psi_0_eq;
T3_store_a = mean(T(T3_node,:)) * Dim_T * psi_0_eq;

T4_store_a  = movmean(T4_store_a, 30/dim_t/(dt*gap));
T3_store_a  = movmean(T3_store_a, 30/dim_t/(dt*gap));
T34_store_a = movmean(T34_store_a, 30/dim_t/(dt*gap));

T4_store_a = T4_store_a(1:30/dim_t/(dt*gap):end);
T3_store_a = T3_store_a(1:30/dim_t/(dt*gap):end);
T34_store_a = T34_store_a(1:30/dim_t/(dt*gap):end);

T4_store_a = T4_store_a(1:30/dim_t/(dt*gap):end);
T3_store_a = T3_store_a(1:30/dim_t/(dt*gap):end);
T34_store_a = T34_store_a(1:30/dim_t/(dt*gap):end);


Starting_year = 325;
Ending_year = 375;
Starting_year = 120;
Ending_year = 170;
dt = 8/24/dim_t/10; % time step, every dt*dim_t hours

LL2 = floor(Starting_year*360/dim_t/(dt*gap)); % Starting time of plotting
LL3 = floor(Ending_year*360/dim_t/(dt*gap)); % Final time of plotting
ll = length(dt*LL2:dt: dt*LL3);
idx = LL2:100:LL2+ll-1;

LL0 = Starting_year*360;
LL1 = Ending_year*360+1;
ll2 = length(LL0:LL1);
idx2 = (LL0:LL1);


[xx,yy] = meshgrid(x*dim_x,idx/360);
%[xx,yy] = meshgrid(120:(280-120)/No,idx/360);

temp_x = (-32:No-1) * dx; % grid points in the x axis
temp_lon = (-32+1:No)*360/Na+120;

[temp_xx,temp_yy] = meshgrid(temp_x*dim_x/1000,idx/360);
[temp_xx_lon,temp_yy] = meshgrid(temp_lon,idx/360);

temp_MJO = MJO([Na-32+1:Na,1:No],:);
temp_A_phy = A_phy([Na-32+1:Na,1:No],:);
temp_u_phy = u_phy([Na-32+1:Na,1:No],:);
temp_theta_phy = theta_phy([Na-32+1:Na,1:No],:);
temp_Q_phy = Q_phy([Na-32+1:Na,1:No],:);
temp_u_bar = u_bar([Na-32+1:Na,1:No],:);
temp_A_bar = A_bar([Na-32+1:Na,1:No],:);
temp_noise = noise_record([Na-32+1:Na,1:No],:);


% Calculate the magnitude of the MJO projection
MJO_magnitude = abs(temp_MJO);
% Define the number of time steps in one year (assuming daily resolution)
steps_per_year = 360;
% Calculate the 1-year running mean of the MJO magnitude
MJO_1yr_mean = movmean(MJO_magnitude, steps_per_year, 2);

pos_x = [0.025:0.09:0.925];
pos_y = 0.18;
pos_w = 0.06;
pos_h = 0.7;



fig = figure('Position',[130.6,152.2,1778.4,791.2]);
colormap jet
ax = subplot('Position',[pos_x(1), pos_y, pos_w, pos_h]);
hold(ax,'on')
contourf(ax, temp_xx_lon, temp_yy, MJO_1yr_mean(:,idx)', 40, 'LineStyle','none')
plot(ax, [120 120], [idx(1) idx(end)]/360, 'k--', 'LineWidth', 1.5)
title(ax,['$(a) |e_{MJO}|$'],'Interpreter','latex')
set(ax,'FontSize',14)
xlabel(ax,'Longitude')
c = colorbar(ax,'southoutside');
set(c,'Position',[pos_x(1) 0.075 pos_w 0.027])
clim(ax,[0 0.4])
set(ax,'Layer','top','Box','on','LineWidth',0.5)


ax = subplot('Position',[pos_x(2),pos_y,pos_w,pos_h]);
hold on
contourf(ax, temp_xx_lon,temp_yy,(temp_A_phy(:,idx))',40,'linestyle','none')
plot(ax, [120,120],[idx(1),idx(end)]/360,'k--','LineWidth',1.5)
title('$(b) a^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim(ax,[-20,20]);
set(c,'position',[pos_x(2),0.075,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(3),pos_y,pos_w,pos_h]);
hold on
%levels=-20:.5:20;
contourf(ax, temp_xx_lon,temp_yy,temp_u_phy(:,idx)',40,'linestyle','none')
plot(ax, [120,120],[idx(1),idx(end)]/360,'k--','LineWidth',1.5)
title('$(c) u^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim(ax,[-20,20]);
set(c,'position',[pos_x(3),0.075,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(4),pos_y,pos_w,pos_h]);
hold on
%contourf(temp_xx,temp_yy, temp_noise(:,idx)',40,'linestyle','none')
contourf(ax, temp_xx_lon,temp_yy,temp_Q_phy(:,idx)',40,'linestyle','none')
title('$(d) q^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
plot(ax, [120,120],ylim,'k--','LineWidth',1.5)
xlabel('Longitude')
c = colorbar('southoutside');
clim([-10,10]);
set(c,'position',[pos_x(4),0.075,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)


subplot('Position',[pos_x(5),pos_y,pos_w,pos_h])
%contourf(temp_xx,temp_yy, (tanh(2*temp_Eq(:,idx).^2).*temp_eta')',40,'linestyle','none')
contourf(temp_xx_lon(:,end-No+1:end),temp_yy(:,end-No+1:end),A_bar(1:No,idx)',40,'linestyle','none')
title('$(e) \bar{a}$','Interpreter','latex')
set(gca,'fontsize',14)
%plot([0,0],[idx(1),idx(end)]/360,'k','LineWidth',1.5)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
%clim([0,1.8]);
set(c,'position',[pos_x(5),0.075,pos_w,0.027])

subplot('Position',[pos_x(6),pos_y,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy(:,end-No+1:end),u_bar(1:No,idx)',40,'linestyle','none')
title('$(f) \bar{u}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
%clim([-3,3]);
set(c,'position',[pos_x(6),0.075,pos_w,0.027])


subplot('Position',[pos_x(7),pos_y,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy(:,end-No+1:end),Hov_U_a(1:No,idx)',40,'linestyle','none')
title('$(g) U$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
set(c,'position',[pos_x(7),0.075,pos_w,0.027])

subplot('Position',[pos_x(8),pos_y,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy(:,end-No+1:end),Hov_H_a(1:No,idx)',40,'linestyle','none')
title('$(h) H$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
set(c,'position',[pos_x(8),0.075,pos_w,0.027])

subplot('Position',[pos_x(9),pos_y,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy(:,end-No+1:end),Hov_T_a(1:No,idx)',40,'linestyle','none')
title('$(i) T$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([-3.5,3.5]);
set(c,'position',[pos_x(9),0.075,pos_w,0.027])


%u_phy2 = eta_temp'.*ap_record*psi_0_eq * dim_u;

temp_u = u_phy+ u_bar;
subplot('Position',[pos_x(10),pos_y,pos_w,pos_h])
plot(mean(temp_u(1:No/2,(LL2:1:LL2+ll-1)),1),(LL2:1:LL2+ll-1)/360)
%plot(ap_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
title('$(j) u_W$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('m·s^{-1}')
ylim([LL2/360,(LL2+ll-1)/360])
%c = colorbar('southoutside');
xlim([-20,20])
%set(c,'position',[pos_x(10),0.075,pos_w,0.027])

subplot('Position',[pos_x(11),pos_y,pos_w,pos_h])
plot(I(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
%plot(qq_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
title('$(k) I$','Interpreter','latex')
set(gca,'fontsize',14)
xlabel(' ')
ylim([LL2/360,(LL2+ll-1)/360])
%c = colorbar('southoutside');
%xlim([-3,3])
%set(c,'position',[pos_x(11),0.075,pos_w,0.027])

%% Figures
fig = figure('Position',[130.6,152.2,1778.4,850]);

annotation('textbox',[0.465012595591544,0.951741996233517,0.147223121907332,0.051789077212806],'String','EP Event','Fontsize',16,'LineStyle','none')
annotation('textbox',[0.465687359424202,0.479519774011298,0.094816464237516,0.051789077212806],'String','CP Event','Fontsize',16,'LineStyle','none')

Starting_year1 = 206.5;
Ending_year1 = 210.5;
LL2_1 = floor((Starting_year1*360+1)/dim_t/(dt*gap));
LL3_1 = floor(Ending_year1*360/dim_t/(dt*gap));
ll_1 = length(dt*LL2_1:dt: dt*LL3_1);
idx_1 = LL2_1:10:LL2_1+ll_1-1;
LL0_1 = Starting_year1*360+1;
LL1_1 = Ending_year1*360;
ll2_1 = length(LL0_1:LL1_1);
idx2_1 = (LL0_1:LL1_1);

Starting_year2 = 1783;
Ending_year2 = 1787;
LL2_2 = floor((Starting_year2*360+1)/dim_t/(dt*gap));
LL3_2 = floor(Ending_year2*360/dim_t/(dt*gap));
ll_2 = length(dt*LL2_2:dt: dt*LL3_2);
idx_2 = LL2_2:10:LL2_2+ll_2-1;
LL0_2 = Starting_year2*360+1;
LL1_2 = Ending_year2*360;
ll2_2 = length(LL0_2:LL1_2);
idx2_2 = (LL0_2:LL1_2);

[temp_xx_lon, temp_yy_1] = meshgrid(temp_lon, idx_1/360);
[~, temp_yy_2] = meshgrid(temp_lon, idx_2/360);


pos_x = [0.03:0.09:0.93];
pos_h = 0.32;
pos_y1 = 0.61;
pos_y2 = 0.16;
pos_w = 0.06;

colormap jet

% ========== 绘制第一个时间段（上面一行）==========
ax = subplot('Position', [pos_x(1), pos_y1, pos_w, pos_h]);
hold(ax, 'on')
contourf(ax, temp_xx_lon, temp_yy_1, temp_MJO(:, idx_1)', 40, 'LineStyle', 'none')
plot(ax, [120 120], [idx_1(1) idx_1(end)]/360, 'k--', 'LineWidth', 1.5)
set(gca,'fontsize',14)
xlabel('Longitude')
title(ax, '$(a) MJO$', 'Interpreter', 'latex')
%c = colorbar(ax, 'southoutside');
%set(c, 'Position', [pos_x(1) pos_y1-0.105 pos_w 0.02])
clim(ax, [-0.7 0.7])
set(ax, 'Layer', 'top', 'Box', 'on', 'LineWidth', 0.5)

ax = subplot('Position',[pos_x(2),pos_y1,pos_w,pos_h]);
hold on
contourf(ax, temp_xx_lon,temp_yy_1,(temp_A_phy(:,idx_1))',40,'linestyle','none')
plot(ax, [120,120],[idx_1(1),idx_1(end)]/360,'k--','LineWidth',1.5)
title('$(b) a^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim(ax,[-20,20]);
%set(c,'position',[pos_x(2),pos_y1-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(3),pos_y1,pos_w,pos_h]);
hold on
%levels=-20:.5:20;
contourf(ax, temp_xx_lon,temp_yy_1,temp_u_phy(:,idx_1)',40,'linestyle','none')
plot(ax, [120,120],[idx_1(1),idx_1(end)]/360,'k--','LineWidth',1.5)
title('$(c) u^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim(ax,[-20,20]);
%set(c,'position',[pos_x(3),pos_y1-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(4),pos_y1,pos_w,pos_h]);
hold on
%contourf(temp_xx,temp_yy, temp_noise(:,idx)',40,'linestyle','none')
contourf(ax, temp_xx_lon,temp_yy_1,temp_Q_phy(:,idx_1)',40,'linestyle','none')
title('$(d) q^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
plot(ax, [120,120],ylim,'k--','LineWidth',1.5)
xlabel('Longitude')
%c = colorbar('southoutside');
clim([-10,10]);
%set(c,'position',[pos_x(4),pos_y1-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)


subplot('Position',[pos_x(5),pos_y1,pos_w,pos_h])
%contourf(temp_xx,temp_yy, (tanh(2*temp_Eq(:,idx).^2).*temp_eta')',40,'linestyle','none')
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_1(:,end-No+1:end),A_bar(1:No,idx_1)',40,'linestyle','none')
title('$(e) \bar{a}$','Interpreter','latex')
set(gca,'fontsize',14)
%plot([0,0],[idx(1),idx(end)]/360,'k','LineWidth',1.5)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim([0,2]);
%set(c,'position',[pos_x(5),pos_y1-0.105,pos_w,0.027])

subplot('Position',[pos_x(6),pos_y1,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_1(:,end-No+1:end),u_bar(1:No,idx_1)',40,'linestyle','none')
title('$(f) \bar{u}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim([-1,1]);
%set(c,'position',[pos_x(6),pos_y1-0.105,pos_w,0.027])


subplot('Position',[pos_x(7),pos_y1,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_1(:,end-No+1:end),Hov_U_a(1:No,idx_1)',40,'linestyle','none')
title('$(g) U$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim([-1.5,1.5]);
%set(c,'position',[pos_x(7),pos_y1-0.105,pos_w,0.027])

subplot('Position',[pos_x(8),pos_y1,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_1(:,end-No+1:end),Hov_H_a(1:No,idx_1)',40,'linestyle','none')
title('$(h) H$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
clim([-50,50]);
%c = colorbar('southoutside');
%set(c,'position',[pos_x(8),pos_y1-0.105,pos_w,0.027])

subplot('Position',[pos_x(9),pos_y1,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_1(:,end-No+1:end),Hov_T_a(1:No,idx_1)',40,'linestyle','none')
title('$(i) T$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
%c = colorbar('southoutside');
clim([-3.5,3.5]);
%set(c,'position',[pos_x(9),pos_y1-0.105,pos_w,0.027])


%u_phy2 = eta_temp'.*ap_record*psi_0_eq * dim_u;

temp_u = u_phy+ u_bar;
subplot('Position',[pos_x(10),pos_y1,pos_w,pos_h])
plot(mean(temp_u(1:No/2,(LL2_1:1:LL2_1+ll_1-1)),1),(LL2_1:1:LL2_1+ll_1-1)/360)
%plot(ap_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
title('$(j) u_W$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('m·s^{-1}')
ylim([LL2_1/360,(LL2_1+ll_1-1)/360])
%c = colorbar('southoutside');
xlim([-15,15])
%set(c,'position',[pos_x(10),0.075,pos_w,0.027])

subplot('Position',[pos_x(11),pos_y1,pos_w,pos_h])
plot(I(LL2_1:1:LL2_1+ll_1-1),(LL2_1:1:LL2_1+ll_1-1)/360)
%plot(qq_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
title('$(k) I$','Interpreter','latex')
set(gca,'fontsize',14)
xlabel(' ')
ylim([LL2_1/360,(LL2_1+ll_1-1)/360])
xlim([0,1])

%%

ax = subplot('Position', [pos_x(1), pos_y2, pos_w, pos_h]);
hold(ax, 'on')
contourf(ax, temp_xx_lon, temp_yy_2, temp_MJO(:, idx_2)', 40, 'LineStyle', 'none')
plot(ax, [120 120], [idx_2(1) idx_2(end)]/360, 'k--', 'LineWidth', 1.5)
%title(ax, '$(l) |e_{MJO}|$', 'Interpreter', 'latex')
set(gca,'fontsize',14)
%ylabel(ax, 'Year', 'FontSize', 12)
xlabel(ax, 'Longitude')
c = colorbar(ax, 'southoutside');
set(c, 'Position', [pos_x(1) pos_y2-0.105 pos_w 0.027])
clim(ax, [-0.7 0.7])
set(ax, 'Layer', 'top', 'Box', 'on', 'LineWidth', 0.5)

ax = subplot('Position',[pos_x(2),pos_y2,pos_w,pos_h]);
hold on
contourf(ax, temp_xx_lon,temp_yy_2,(temp_A_phy(:,idx_2))',40,'linestyle','none')
plot(ax, [120,120],[idx_2(1),idx_2(end)]/360,'k--','LineWidth',1.5)
%title('$(b) a^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim(ax,[-20,20]);
set(c,'position',[pos_x(2),pos_y2-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(3),pos_y2,pos_w,pos_h]);
hold on
%levels=-20:.5:20;
contourf(ax, temp_xx_lon,temp_yy_2,temp_u_phy(:,idx_2)',40,'linestyle','none')
plot(ax, [120,120],[idx_2(1),idx_2(end)]/360,'k--','LineWidth',1.5)
%title('$(c) u^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim(ax,[-20,20]);
set(c,'position',[pos_x(3),pos_y2-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

ax = subplot('Position',[pos_x(4),pos_y2,pos_w,pos_h]);
hold on
%contourf(temp_xx,temp_yy, temp_noise(:,idx)',40,'linestyle','none')
contourf(ax, temp_xx_lon,temp_yy_2,temp_Q_phy(:,idx_2)',40,'linestyle','none')
%title('$(d) q^{\prime}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
plot(ax, [120,120],ylim,'k--','LineWidth',1.5)
xlabel('Longitude')
c = colorbar('southoutside');
clim([-10,10]);
set(c,'position',[pos_x(4),pos_y2-0.105,pos_w,0.027])
set(ax,'Layer','top','Box','on','LineWidth',0.5)


subplot('Position',[pos_x(5),pos_y2,pos_w,pos_h])
%contourf(temp_xx,temp_yy, (tanh(2*temp_Eq(:,idx).^2).*temp_eta')',40,'linestyle','none')
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_2(:,end-No+1:end),A_bar(1:No,idx_2)',40,'linestyle','none')
%title('$(e) \bar{a}$','Interpreter','latex')
set(gca,'fontsize',14)
%plot([0,0],[idx(1),idx(end)]/360,'k','LineWidth',1.5)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([0,2]);
set(c,'position',[pos_x(5),pos_y2-0.105,pos_w,0.027])

subplot('Position',[pos_x(6),pos_y2,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_2(:,end-No+1:end),u_bar(1:No,idx_2)',40,'linestyle','none')
%title('$(f) \bar{u}$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([-1,1]);
set(c,'position',[pos_x(6),pos_y2-0.105,pos_w,0.027])


subplot('Position',[pos_x(7),pos_y2,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_2(:,end-No+1:end),Hov_U_a(1:No,idx_2)',40,'linestyle','none')
%title('$(g) U$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([-1.5,1.5]);
set(c,'position',[pos_x(7),pos_y2-0.105,pos_w,0.027])

subplot('Position',[pos_x(8),pos_y2,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_2(:,end-No+1:end),Hov_H_a(1:No,idx_2)',40,'linestyle','none')
%title('$(h) H$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([-50,50]);
set(c,'position',[pos_x(8),pos_y2-0.105,pos_w,0.027])

subplot('Position',[pos_x(9),pos_y2,pos_w,pos_h])
contourf(temp_xx_lon(:,end-No+1:end),temp_yy_2(:,end-No+1:end),Hov_T_a(1:No,idx_2)',40,'linestyle','none')
%title('$(i) T$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('Longitude')
c = colorbar('southoutside');
clim([-3.5,3.5]);
set(c,'position',[pos_x(9),pos_y2-0.105,pos_w,0.027])


%u_phy2 = eta_temp'.*ap_record*psi_0_eq * dim_u;

temp_u = u_phy+ u_bar;
subplot('Position',[pos_x(10),pos_y2,pos_w,pos_h])
plot(mean(temp_u(1:No/2,(LL2_2:1:LL2_2+ll_2-1)),1),(LL2_2:1:LL2_2+ll_2-1)/360)
%plot(ap_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
%title('$(j) u_W$','Interpreter','latex')
set(gca,'fontsize',14)
%ylim([300,Time_total])
xlabel('m·s^{-1}')
ylim([LL2_2/360,(LL2_2+ll_2-1)/360])
%c = colorbar('southoutside');
xlim([-15,15])
%set(c,'position',[pos_x(10),0.105,pos_w,0.027])

subplot('Position',[pos_x(11),pos_y2,pos_w,pos_h])
plot(I(LL2_2:1:LL2_2+ll_2-1),(LL2_2:1:LL2_2+ll_2-1)/360)
%plot(qq_record(LL2:1:LL2+ll-1),(LL2:1:LL2+ll-1)/360)
%title('$(k) I$','Interpreter','latex')
set(gca,'fontsize',14)
xlabel(' ')
ylim([LL2_2/360,(LL2_2+ll_2-1)/360])
xlim([0,1])
