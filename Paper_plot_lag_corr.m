Starting_year = 100;
Ending_year = 2000;

ifte = 2;

LL = 10;
lag = -3*360/LL:360/LL*2; % Lag in months

dim_t = 34;
dt = 8/24/dim_t/10; % time step, every dt*dim_t hours

gap = 30;

LL2 = floor((Starting_year*360+1)/dim_t/(dt*gap)); % Starting time of plotting
LL3 = floor(Ending_year*360/dim_t/(dt*gap)); % Final time of plotting
ll = length(dt*LL2:dt: dt*LL3);
idx = LL2:LL:LL2+ll-1;

idx2 = idx;


idx = idx2(1)-lag(1)*LL:LL:idx2(end)-lag(end)*LL;


west_node = 1:round(No*1/2);
east_node = round(No*1/2)+1:No;
T3_node = 31:53;
T4_node = 14:31;
T34_node = 25:43;

T4_store_a = mean(T(T4_node,:)) * Dim_T * psi_0_eq;
Tw_store_a = mean(T(west_node,:)) * Dim_T * psi_0_eq;
T34_store_a = mean(T(T34_node,:)) * Dim_T * psi_0_eq;
Te_store_a = mean(T(east_node,:)) * Dim_T * psi_0_eq;
Tc_store_a = mean(T(T4_node,:)) * Dim_T * psi_0_eq;
T3_store_a = mean(T(T3_node,:)) * Dim_T * psi_0_eq;


temp_x = (-32:No-1) * dx; % grid points in the x axis
temp_x2 = 120-160/55*32:160/55:280;

[temp_xx,temp_yy] = meshgrid(temp_x*dim_x/1000,idx2/360);
[temp_xx2,temp_yy2] = meshgrid(temp_x2,idx2/360);

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


if ifte == 0
    T_E_lag = Tc_store_a(idx)';
elseif ifte == 1
    T_E_lag = Te_store_a(idx)';
elseif ifte == 2
    T_E_lag = T34_store_a(idx)';
end
MJO_lag = temp_MJO(:,idx2)';
abs_MJO_lag = abs(MJO_lag);
u_lag = temp_u_phy(:,idx2)';
abs_u_lag = abs(u_lag);
u_bar_lag = u_bar(1:No,idx2)';
H_lag = Hov_H_a(:,idx2)';
U_lag = Hov_U_a(:,idx2)';
T_lag = Hov_T_a(:,idx2)';
a_lag = temp_A_phy(:,idx2)';
a_bar_lag = A_bar(1:No,idx2)';



reg_e_MJO = lagged_corr(MJO_lag, T_E_lag, lag);
reg_abs_e_MJO = lagged_corr(abs_MJO_lag, T_E_lag, lag);
reg_u_prime = lagged_corr(u_lag, T_E_lag, lag);
reg_abs_u_prime = lagged_corr(abs_u_lag, T_E_lag, lag);
reg_u_bar = lagged_corr(u_bar_lag, T_E_lag, lag);
reg_H = lagged_corr(H_lag, T_E_lag, lag);
reg_U = lagged_corr(U_lag, T_E_lag, lag);
reg_T = lagged_corr(T_lag, T_E_lag, lag);
reg_a = lagged_corr(a_lag, T_E_lag, lag);
reg_a_bar = lagged_corr(a_bar_lag, T_E_lag, lag);

% Create the figure
figure
colormap jet
subplot(1,9,1)
contourf(temp_x2, lag/360*LL,reg_e_MJO,'LineStyle','none')
title('(a) $MJO$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,2)
contourf(temp_x2, lag/360*LL, reg_abs_e_MJO,'LineStyle','none')
title('(b) $|MJO|$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,3)
contourf(temp_x2, lag/360*LL, reg_a,'LineStyle','none')
title('(c) $a^{\prime}$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,4)
contourf(temp_x2, lag/360*LL, reg_u_prime,'LineStyle','none')
title('(d) $u^{\prime}$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,5)
contourf(temp_x2, lag/360*LL, reg_abs_u_prime,'LineStyle','none')
title('(e) $|u^{\prime}|$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,6)
contourf(120:160/55:280, lag/360*LL, reg_u_bar,'LineStyle','none')
title('(f) $\bar{u}$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14, 'YDir', 'normal')
colorbar('southoutside')



subplot(1,9,7)
contourf(120:160/55:280, lag/360*LL, reg_U,'LineStyle','none')
title('(g) $U$','Interpreter','latex')
ylabel('Lag (months)')
xlabel('Longitude')
set(gca, 'FontSize',14, 'YDir', 'normal')
colorbar('southoutside')

subplot(1,9,8)
contourf(120:160/55:280, lag/360*LL, reg_H,'LineStyle','none')
title('(h) $H$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

subplot(1,9,9)
contourf(120:160/55:280, lag/360*LL, reg_T,'LineStyle','none')
title('(i) $T$','Interpreter','latex')
ylabel('Lag (years)')
xlabel('Longitude')
set(gca, 'FontSize',14)
colorbar('southoutside')

% Add a red line at x = 0
for i = 1:5
    subplot(1,9,i)
    hold on
    plot([120 120], ylim, 'r-', 'LineWidth', 2)
    hold off
end

for i = 1:9
    subplot(1,9,i)
    hold on
    plot(xlim, [0,0], 'k--', 'LineWidth', 2)
    hold off
end

function corr_matrix = lagged_corr(X, y, lag)
    [~, n_x] = size(X);
    corr_matrix = zeros(length(lag), n_x);
    for i = 1:length(lag)
        X_lag = X(-lag(1)+lag(i)+1:-lag(1)+lag(i)+length(y), :);
        y_lag = y;
        for j = 1:n_x
            corr_matrix(i,j) = corr(X_lag(:,j), y_lag, 'rows', 'complete');
        end
    end
end