load uwnd_new_data.mat
load hgt_new_data.mat
load Obs_fine_2020.mat
load obs_MJO_new.mat  % Load MJO data separately

%% Setup parameters and coordinate systems
Left1_end = (120) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST
Middle1_end = 200 / 2.5 * 10; % Middle of Pacific ocean 160W for SST
Nino1_end = 190 / 2.5 * 10; % Middle of Pacific ocean 170W for SST
Nino2_end = 240 / 2.5 * 10; % Middle of Pacific ocean 120W for SST

Left3_end = (120) / 2.5; % Pacific ocean left boundary 120E for wind bursts
Right3_end = 280 / 2.5; % Pacific ocean right boundary 80W for wind bursts
Middle3_end = 200 / 2.5; % Middle of Pacific ocean 160W for wind bursts
Left3_1 = (25+2.5)/2.5;
Left1_C = 160/2.5*10;
Right1_C = 210/2.5*10;

% Calculate physical wind from modes
psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator

K_total = (uwnd_mode_0_rmmean_3modes - hgt_mode_0_rmmean_3modes)/2;
R_total = -(uwnd_mode_0_rmmean_3modes + hgt_mode_0_rmmean_3modes)/4 + (uwnd_mode_2_rmmean_3modes - hgt_mode_2_rmmean_3modes)/2/sqrt(2);
K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/2;
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/4 + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes)/2/sqrt(2);

u_total_obs = ( (K_total - R_total) * psi_0 + 1/sqrt(2) * R_total * psi_2 ) * 50; % unit 50m/s
u_phy_obs = ( (K_3modes - R_3modes) * psi_0 + 1/sqrt(2) * R_3modes * psi_2 ) * 50; % unit 50m/s
u_bar_obs = u_total_obs - u_phy_obs;

load bigSST_oneyear_obs.mat
load bigu_oneyear_obs.mat
SST_b = repmat(T_SC_obs,20,1)';
u_b = repmat(u_SC_obs,20,1)';

% Extract 20 years of data
years_to_extract = 20;
days_per_year = 365;
start_year = 1982;
total_days = years_to_extract * days_per_year;
time1 = (start_year-1982)*365+1:(start_year-1982)*365+total_days;
time2 = (start_year-1982)*12+1:(start_year-1982)*12+total_days/30;
time3 = time1+1095;
ttt= Y;

MJO_temp = MJO(:, time3);
wind = u_phy_obs(:, time3);
SST = sst_a_fine(time1, Left1_end:Right1_end)' + SST_b;

% Calculate Nino 3.4 equivalent (Eastern Pacific SST anomaly)
nino34_indices = Nino1_end:Nino2_end;
T34_new = mean(sst_a_fine(time1, nino34_indices), 2)';
T34_new = movmean(T34_new, 90); % 90-day running mean
T34_new = T34_new(1:30:end);

Event_count = T34_new > 0.5;
diff_array = diff(Event_count);
Eve_start = find(diff_array == 1);  % ENSO event start
Eve_end = find(diff_array == -1);   % ENSO event end

% If the last event hasn't ended, set to end of time series
if length(Eve_start) > length(Eve_end)
    Eve_end = [Eve_end, length(Event_count)];
end

% Initialize MJO classification
MJO_count = zeros(size(Event_count));
MJO_count(Event_count == 1) = 1;  % 1 = during El Niño

% Mark MJO events 9 months prior to each ENSO event
for i = 1:length(Eve_start)
    % Define time window for 9 months before the event
    pre_event_start = max(1, Eve_start(i) - 9);
    pre_event_end = Eve_start(i) - 1;
    
    % Mark time points within 9 months before the event
    if pre_event_end >= pre_event_start
        MJO_count(pre_event_start:pre_event_end) = 2;  % 2 = prior to El Niño
    end
end

% Handle overlap cases: if a time point is both during ENSO event and prior period,
% prioritize marking as during event
MJO_count(Event_count == 1) = 1;

%% figure1
MJO_x = size(MJO,1);

[xx_new,yy_new] = meshgrid(X_sst(Left1_end:Right1_end),ttt(time1));
[xx_new2,yy_new2] = meshgrid(40:161/56:280,time1/365);
[xx_new3,yy_new3] = meshgrid(X_uwnd(Left3_end:Right3_end),ttt(time1));
temp_MJO2 = MJO(Left3_1:Right3_end,time3);
% Create the figure
figure('Position', [593.8,7.4,1384,1032.8]);
% Subplot (a): Warm Pool Edge
subplot('Position',[0.05,0.15,0.14,0.815])
hold on
%contourf(X_uwnd(Left3_1:Right3_end), lag/365*LL, reg_MJO,'LineStyle','none')
[CC1,hh1] = contour(X_uwnd(Left3_1:Right3_end),ttt(time1),temp_MJO2', [0.06 0.06], 'LineColor', 'r', 'LineWidth', 1.5,'Visible','off');

%[CC1,hh1] = contour(xx_new2,yy_new2,temp_MJO2', [0.2 0.2], 'LineColor', 'r', 'LineWidth', 1.5,'Visible','off');
%contour(xx_new,yy_new,MJO_temp(1:No,:)', [-0.19 -0.19], 'LineColor', 'b', 'LineWidth', 1.5);
%contourf(xx_new,yy_new,MJO_temp(1:No,:)','LineStyle','none')

    % [xx,yy] = meshgrid(X_sst(Left1_end:Right1_end),ttt(time2));
    % contourf(xx,yy,sst_a_fine(time2,Left1_end:Right1_end),'LineStyle','none')

[wpee1, h1] = contour(X_sst(Left1_end:Right1_end),ttt(time1), SST', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5,'Visible','off');
[xa1, ya1] = extractLargestContour(wpee1);
plot(xa1, ya1,'m', 'LineWidth', 1.5)
%[wpee1, h1] = contour(xx_new(:,6:end-15),yy_new(:,6:end-15), SST(6:end-15,:)', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5);
[wpee2, h2] = contour(xx_new(:,60:end),yy_new(:,60:end), SST_b(60:end,:)', [28.5 28.5], 'LineColor', 'b', 'LineWidth', 1.5,'Linestyle','--');
box on
xlabel('Longitude')
grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
title('(a) MJO events','FontWeight','normal')
ylabel('years');
set(gca,'FontSize',14)
contourData = extractContourData(CC1);
clear CC1 hh1
for i = 1:length(contourData)
    level = contourData(i).level;
    x = contourData(i).x;
    y = contourData(i).y;

    time_MJO = mean(y);
    kind_MJO = MJO_count(ceil(time_MJO*12-start_year*12));

    %MJO_color = ['b','r','m'];
    MJO_color = [0.5,0.5,0.5;1,0,0;0.93,0.69,0.13];

    fprintf('Contour level: %.2f, Number of points: %d\n', level, length(x));

    plot(x, y,'color',[MJO_color(kind_MJO+1,:)],'MarkerFaceColor',[MJO_color(kind_MJO+1,:)],'LineWidth',1.5);
end

[xa2, ya2] = extractLargestContour(wpee2);
%clear wpee1 h1 wpee2 h2

[~, ia, ~] = unique(ya1, 'first');
ya1 = ya1(ia);
xa1 = xa1(ia);
[~, ia2, ~] = unique(ya2, 'first');
ya2 = ya2(ia2);
xa2 = xa2(ia2);

ttt2 = ttt(time1);
dist_new = linspace(ttt2(1), ttt2(end), total_days);
xa1_new= interp1(ya1, xa1, dist_new, 'linear');
xa2_new= interp1(ya2, xa2, dist_new, 'linear');

%%
subplot('Position',[0.25,0.15,0.14,0.815])
data = xa1_new-xa2_new;
positive_data = max(data, 0);
negative_data = min(data, 0);

% Create area plot
hold on;
area(ttt(time1), positive_data*111/1000, 'FaceColor', [1 0.7 1], 'EdgeColor', 'k', 'LineWidth', 1);  % Light pink
area(ttt(time1), negative_data*111/1000, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'k', 'LineWidth', 1);  % Light blue
plot(ttt(time1), data*111/1000, 'k-', 'LineWidth', 1);

ylim([-6, 9]);
ylabel('10^3 km');

title('(b) Warm pool anomaly','FontWeight','normal');

grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');

set(gca, 'Box', 'on');
set(gca, 'XAxisLocation', 'bottom');  % Move X-axis to bottom
set(gca,'FontSize',14)

%plot([0 0], ylim, 'k-', 'LineWidth', 1);

view([90 -90]);

LL = 1;
%%
ax = subplot('Position',[0.45,0.15,0.14,0.815]);
hold on
colormap jet
contourf(ax,xx_new3(1:LL:end,:),yy_new3(1:LL:end,:),(u_total_obs(Left3_end:Right3_end,time3(1:LL:end)))'+u_b(:,1:LL:end)','LineStyle','none')
%contour(xx_new(:,6:end),yy_new(:,6:end), SST(6:end,:)', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5);
[xa1, ya1] = extractLargestContour(wpee1);
plot(xa1, ya1,'m', 'LineWidth', 1.5)
contour(xx_new(:,60:end),yy_new(:,60:end), SST_b(60:end,:)', [28.5 28.5], 'LineColor', 'b', 'LineWidth', 1.5,'Linestyle','--');
box on
xlim([120,240])
set(gca,'xtick',120:60:240);
xlabel('Longitude')
title('(c) Wind','FontWeight','normal');
set(gca,'FontSize',14)
colorbar('southoutside','Position',[0.450289017341041,0.063774861617941,0.139884393063584,0.020655822359928])
clim([-25,20])
set(gca,'linewidth',1);
ylim([ttt2(1),ttt2(end)])
set(ax,'Layer','top','Box','on','LineWidth',0.5)

%% figure4
subplot('Position',[0.65,0.15,0.14,0.815])
dataset = T34_new;
flag1 = 1;
flag2 = 1;
days2 = ttt2(1:30:end);
xn = days2(dataset < 0);
yn = dataset(dataset < 0);
xp = days2(dataset >= 0);
yp = dataset(dataset >= 0);

light_red = [255, 182, 193] / 255;
dark_red = [220, 20, 60] / 255;
light_blue = [173, 216, 230] / 255;
dark_blue = [0, 0, 139] / 255;

hold on

for i = 1:length(dataset)
    if dataset(i) >= 0
        if dataset(i) <= 0.5
            plot([0, dataset(i)], [days2(i), days2(i)], 'Color', light_red, 'LineWidth', 2);
        else
            plot([0, 0.5], [days2(i), days2(i)], 'Color', light_red, 'LineWidth', 2);
            plot([0.5, dataset(i)], [days2(i), days2(i)], 'Color', dark_red, 'LineWidth', 2);
        end
    else
        if dataset(i) >= -0.5
            plot([dataset(i), 0], [days2(i), days2(i)], 'Color', light_blue, 'LineWidth', 2);
        else
            plot([-0.5, 0], [days2(i), days2(i)], 'Color', light_blue, 'LineWidth', 2);
            plot([dataset(i), -0.5], [days2(i), days2(i)], 'Color', dark_blue, 'LineWidth', 2);
        end
    end
end

% Set grid
grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
box on
ylim([ttt2(1),ttt2(end)])
title('(d) Nino 3.4','FontWeight','normal')
set(gca,'FontSize',14)
%%
subplot('Position',[0.85,0.15,0.14,0.815])
%[xx,yy] = meshgrid(t_model,lon_sst);
contourf(X_sst(Left1_end:Right1_end),ttt(time1),sst_a_fine(time1, Left1_end:Right1_end),30,'linestyle','none')

T_E = mean(sst_a_fine(:,Middle1_end:Right1_end),2);% SST averaged over the eastern Pacific
T_C = mean(sst_a_fine(:,Left1_C:Right1_C),2);% SST averaged over the central Pacific

T_E = movmean(T_E, 90); % 90-day running mean
T_E_3R = T_E(1:30:end)';
T_C = movmean(T_C, 90); % 90-day running mean
T_C_3R = T_C(1:30:end)';

% T_C_3R = T_E';
% T_E_3R = T_C';

caxis([-3,3])
hold on
window = 12;
total_loop = (years_to_extract-1) * 12/window;
colormap(jet)
range_model = time2;
t_model = 1982+range_model/12;
plot([180 180],[t_model(1) t_model(end)],'m--','linewidth',2);
for k = 1:total_loop
    if mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 1.0 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
        plot([120,120],[t_model(1-4+k*window),t_model(1+1+window*k)],'r','linewidth',10)
    elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
        plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'m','linewidth',10)
    elseif mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
        plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'color',[255 97 0]/255,'linewidth',10)
    elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5 || mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5
        plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'b','linewidth',10)
    end
end
%     colorbar
set(gca,'xlim',[120 280]);
set(gca,'xtick',120:60:280);
ylim([ttt2(1),ttt2(end)])

xlabel('longitude');
set(gca,'linewidth',1);
title('(e) SST anomaly','FontWeight','normal');
set(gca,'FontSize',14)
colorbar('southoutside','Position',[0.85028901734104,0.063774861617941,0.139884393063584,0.020655822359928])