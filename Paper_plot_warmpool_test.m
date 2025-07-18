load bigSST_oneyear.mat
load bigu_oneyear.mat
SST_b = repmat(T_SC_one,20,1)';
u_b = repmat(u_SC_new,20,1)';

% Extract 20 years of data
years_to_extract = 20;
days_per_year = 360;
start_year = 724;
%start_year = 124;
total_days = years_to_extract * days_per_year;
time1 = start_year*360+1:start_year*360+total_days;
time2 = start_year*12+1:start_year*12+total_days/30;

MJO_temp = MJO(:, time1);
wind = u_phy(:, time1);
SST = Hov_T_a(:, time1) + SST_b;

T34_new = mean(T(T34_node,time1)) * Dim_T * psi_0_eq;
T34_new = movmean(T34_new, 90);
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

[xx_new,yy_new] = meshgrid(120:161/56:280,time1/360);
[xx_new2,yy_new2] = meshgrid(40:161/56:280,time1/360);
temp_MJO2 = MJO_temp([Na-28+1:Na,1:No],:);
% Create the figure
figure('Position', [593.8,7.4,1384,1032.8]);
% Subplot (a): Warm Pool Edge
subplot('Position',[0.05,0.15,0.14,0.815])
hold on
[CC1,hh1] = contour(xx_new2,yy_new2,temp_MJO2', [0.25 0.25], 'LineColor', 'r', 'LineWidth', 1.5,'Visible','off');
%contour(xx_new,yy_new,MJO_temp(1:No,:)', [-0.19 -0.19], 'LineColor', 'b', 'LineWidth', 1.5);
%contourf(xx_new,yy_new,MJO_temp(1:No,:)','LineStyle','none')
[wpee1, h1] = contour(xx_new(:,6:end),yy_new(:,6:end), SST(6:end,:)', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5,'Visible','off');
[xa1, ya1] = extractLargestContour(wpee1);
plot(xa1, ya1,'m', 'LineWidth', 1.5)
%[wpee1, h1] = contour(xx_new(:,6:end-15),yy_new(:,6:end-15), SST(6:end-15,:)', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5);
[wpee2, h2] = contour(xx_new(:,10:end),yy_new(:,10:end), SST_b(10:end,:)', [28.5 28.5], 'LineColor', 'b', 'LineWidth', 1.5,'Linestyle','--');
box on
xlabel('Longitude')
grid on;
ylim([time1(1)/360,time1(end)/360])
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

dist_new = linspace(time1(1)/360, time1(end)/360, total_days);
xa1_new= interp1(ya1, xa1, dist_new, 'linear');
xa2_new= interp1(ya2, xa2, dist_new, 'linear');

%%
subplot('Position',[0.25,0.15,0.14,0.815])
data = xa1_new-xa2_new;
positive_data = max(data, 0);
negative_data = min(data, 0);

% Create area plot
hold on;
area(time1/360, positive_data*111/1000, 'FaceColor', [1 0.7 1], 'EdgeColor', 'k', 'LineWidth', 1);  % Light pink
area(time1/360, negative_data*111/1000, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'k', 'LineWidth', 1);  % Light blue
plot(time1/360, data*111/1000, 'k-', 'LineWidth', 1);

ylim([-6, 9]);
ylabel('10^3 km');
title('(b) Warm pool anomaly','FontWeight','normal');

grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');

set(gca, 'Box', 'on');
set(gca, 'XAxisLocation', 'bottom');  % Move X-axis to bottom
set(gca,'FontSize',14)

%plot([0 0], ylim, 'k-', 'LineWidth', 1);
xlim([time1(1)/360,time1(end)/360])
view([90 -90]);
LL = 1;
%%
ax = subplot('Position',[0.45,0.15,0.14,0.815]);
hold on
colormap jet
contourf(xx_new(1:LL:end,:),yy_new(1:LL:end,:),(u_phy(1:No,time1(1:LL:end))+u_bar(1:No,time1(1:LL:end)))'+u_b(:,1:LL:end)','LineStyle','none')
%contour(xx_new(:,6:end),yy_new(:,6:end), SST(6:end,:)', [28.5 28.5], 'LineColor', 'k', 'LineWidth', 1.5);
[xa1, ya1] = extractLargestContour(wpee1);
plot(xa1, ya1,'m', 'LineWidth', 1.5)
contour(xx_new(:,10:end),yy_new(:,10:end), SST_b(10:end,:)', [28.5 28.5], 'LineColor', 'b', 'LineWidth', 1.5,'Linestyle','--');
box on
xlim([120,240])
set(gca,'xtick',120:60:240);
xlabel('Longitude')
title('(c) Wind','FontWeight','normal');
set(gca,'FontSize',14)
colorbar('southoutside','Position',[0.450289017341041,0.063774861617941,0.139884393063584,0.020655822359928])
clim([-25,20])
set(gca,'linewidth',1);
ylim([time1(1)/360,time1(end)/360])
set(ax,'Layer','top','Box','on','LineWidth',0.5)
%% figure4
subplot('Position',[0.65,0.15,0.14,0.815])
dataset = T34_new;  % Dataset for plotting
flag1 = 1;
flag2 = 1;
days = time1(1:30:end)/360;
xn = days(dataset < 0);
yn = dataset(dataset < 0);
xp = days(dataset >= 0);
yp = dataset(dataset >= 0);

% Define colors
light_red = [255, 182, 193] / 255;
dark_red = [220, 20, 60] / 255;
light_blue = [173, 216, 230] / 255;
dark_blue = [0, 0, 139] / 255;

hold on

for i = 1:length(dataset)
    if dataset(i) >= 0
        if dataset(i) <= 0.5
            plot([0, dataset(i)], [days(i), days(i)], 'Color', light_red, 'LineWidth', 2);
        else
            plot([0, 0.5], [days(i), days(i)], 'Color', light_red, 'LineWidth', 2);
            plot([0.5, dataset(i)], [days(i), days(i)], 'Color', dark_red, 'LineWidth', 2);
        end
    else
        if dataset(i) >= -0.5
            plot([dataset(i), 0], [days(i), days(i)], 'Color', light_blue, 'LineWidth', 2);
        else
            plot([-0.5, 0], [days(i), days(i)], 'Color', light_blue, 'LineWidth', 2);
            plot([dataset(i), -0.5], [days(i), days(i)], 'Color', dark_blue, 'LineWidth', 2);
        end
    end
end

% Set grid
grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
box on
ylim([time1(1)/360,time1(end)/360])
title('(d) Nino 3.4','FontWeight','normal')
set(gca,'FontSize',14)
%%
subplot('Position',[0.85,0.15,0.14,0.815])
%[xx,yy] = meshgrid(t_model,lon_sst);
contourf(xx_new,yy_new,Hov_T_a(:,time1)',30,'linestyle','none')
T_C_3R = T4_store_a';
T_E_3R = T3_store_a';

caxis([-3,3])
hold on
window = 12;
total_loop = (years_to_extract-1) * 12/window;
colormap(jet)
range_model = time2;
t_model = range_model/12;
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
if kk == 1
    ylabel('model year');
end
ylim([time1(1)/360,time1(end)/360])
xlabel('longitude');
set(gca,'linewidth',1);
title('(e) SST anomaly','FontWeight','normal');
set(gca,'FontSize',14)
colorbar('southoutside','Position',[0.85028901734104,0.063774861617941,0.139884393063584,0.020655822359928])