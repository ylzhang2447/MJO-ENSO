load uwnd_new_data.mat
load hgt_new_data.mat
load Obs_fine_2020.mat
load obs_MJO_new.mat

Left3_end = (120+2.5) / 2.5; % Pacific ocean left boundary 120E for wind bursts
Right3_end = 280 / 2.5; % Pacific ocean right boundary 80W for wind bursts

% Finding the indices of the Pacific ocean domain in observational data sets
Left1_end = (120+2.5) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST
Middle1_end = 200 / 2.5 * 10; % Middle of Pacific ocean 160W for SST

psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator
%Total_WB = -uwnd_mean + uwnd_minus + uwnd_plus; % Total wind bursts time series
K_total = (uwnd_mode_0_rmmean_3modes - hgt_mode_0_rmmean_3modes)/2;
R_total = -(uwnd_mode_0_rmmean_3modes + hgt_mode_0_rmmean_3modes)/4 + (uwnd_mode_2_rmmean_3modes - hgt_mode_2_rmmean_3modes)/2/sqrt(2);

K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/2;
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/4 + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes)/2/sqrt(2);

u_total_obs = ( (K_total - R_total) * psi_0 + 1/sqrt(2) * R_total * psi_2 ) * 50; % unit 50m/s
u_phy_obs = ( (K_3modes - R_3modes) * psi_0 + 1/sqrt(2) * R_3modes * psi_2 ) * 50; % unit 50m/s
u_bar_obs = u_total_obs-u_phy_obs;

time2 = (1982-1982)*365+1:(1987-1982)*365+1; %use these indices for wind, SST, and thermocline
time3 = 1095+time2; %use these indicies for MJO

num = 1;
time_total = 44;
lon = 120+2.5:2.5:280;

%% Wind Event Detect
%% === Step 0: Prepare Intraseasonal Wind Data and Time ===
%u_filt = bandpass(u_phy_obs, [1/100 1/20], 2);  % filter along time dim (dim=2)
%u_filt = u_phy_obs;


% Fourier filtering parameters
time_high_period = 90; % high period cutoff (days) - upper limit of intraseasonal oscillation
time_low_period = 30; % low period cutoff (days) - lower limit of intraseasonal oscillation
max_wavenumber = 10; % maximum spatial wavenumber, keep wavenumbers from 0 to max_wavenumber

dt = 1; % single integration step (day)

total_steps = size(u_phy_obs,2);
% create time axis (in days)
days_per_year = 365; % based on 365-day calendar
total_years = floor(total_steps / days_per_year);
time_axis = linspace(1, total_years, total_steps);

% determine time range to process
t_idx = 1:total_steps;
selected_time = time_axis(t_idx);
time_steps2 = length(t_idx);
time_steps = total_steps;

%fprintf('Original data dimensions: [%d, %d]\n', size(u_phy));
%fprintf('Selected time range: %.2f - %.2f years\n', yr0, yr1);
%fprintf('Processing time steps: %d\n', time_steps);

% extract data for the time range to be processed
u_slice = u_phy_obs;

fprintf('Starting 2D Fourier spatiotemporal filtering...\n');
fprintf('Original data dimensions: [%d, %d]\n', size(u_slice,1), size(u_slice,2));

%% ---------- 1. 2D Fourier Transform -------------------------------------------
% calculate 2D Fourier transform
fprintf('Calculating 2D Fourier transform...\n');

% get actual data dimensions
[n_space, n_time] = size(u_slice);

% calculate 2D FFT
u_fft2d = fft2(u_slice);

% calculate frequency axes
fs_time = 1; % temporal sampling frequency (one point per day)
time_freqs = fs_time * (0:(n_time-1)) / n_time; % temporal frequency axis
spatial_wavenumbers = 0:(n_space-1); % spatial wavenumber axis

%% ---------- 2. Design 2D Filter -------------------------------------------
fprintf('Designing 2D filter...\n');

% temporal frequency filter parameters
time_high_freq = 1/time_high_period; % low frequency corresponding to high period
time_low_freq = 1/time_low_period; % high frequency corresponding to low period

% initialize filter
time_filter = zeros(n_space, n_time);
space_filter = ones(n_space, n_time);

%% Temporal direction filter
for i = 1:n_space
    for j = 1:n_time
        f = time_freqs(j);
        
        % handle frequency periodicity
        if (f >= time_high_freq && f <= time_low_freq) || ...
           (f >= (1-time_low_freq) && f <= (1-time_high_freq))
            time_filter(i, j) = 1;
        end
    end
end

%% Spatial direction filter - using smooth transition
space_transition_width = 2; % transition band width (use smaller transition band for u)

for i = 1:n_space
    for j = 1:n_time
        k = spatial_wavenumbers(i);
        
        % handle low wavenumber region
        if k <= max_wavenumber - space_transition_width
            space_filter(i, j) = 1; % complete pass
        elseif k <= max_wavenumber
            % smooth transition region
            weight = 0.5 * (1 + cos(pi * (k - max_wavenumber + space_transition_width) / space_transition_width));
            space_filter(i, j) = weight;
        
        % handle high wavenumber region (considering FFT periodicity)
        elseif k >= (n_space - max_wavenumber) + space_transition_width
            space_filter(i, j) = 1; % complete pass
        elseif k >= (n_space - max_wavenumber)
            % smooth transition region
            weight = 0.5 * (1 + cos(pi * (n_space - max_wavenumber + space_transition_width - k) / space_transition_width));
            space_filter(i, j) = weight;
        end
        % other regions remain 0 (stop band)
    end
end

% combine temporal and spatial filters
filter_2d = time_filter .* space_filter;

% apply filter
u_fft2d_filtered = u_fft2d .* filter_2d;


% inverse transform
u_filtered = real(ifft2(u_fft2d_filtered));

%u_new = u_phy;
u_filt = u_filtered;

%u_filt = u_phy_obs;

%u_filt = movmean(u_phy_obs, 30, 2, 'omitnan', 'Endpoints', 'fill');
intraseasonal_wind = u_filt(Left3_end:Right3_end,:);  % [time × lon]
ttt = 1979+[1:365*time_total]/365;              % decimal year time vector
[nlon,ntime] = size(intraseasonal_wind);

u_std = std(u_phy_obs(:), 'omitnan');
threshold = 5;  % e.g., 95th percentile filtering
% Event detection parameters
%threshold = 2;            % Wind stress threshold (N/m²)
min_duration = 3;         % Minimum duration (days)
min_zonal_extent = 10;     % Minimum zonal span (degrees)
se = strel('rectangle', [1, 3]);  % Structuring element for morphological closing

% Convert decimal year to datetime
%time_all = datetime(floor(ttt), 1, 1) + calmonths(round(12 * (ttt - floor(ttt))));
time_all = datetime(1979, 1, 1) + days(0:ntime-1);

%% === Step 1: Binary Event Masking and Morphological Closing ===
wwe_regions = intraseasonal_wind > threshold;    % Westerly wind events
ewe_regions = intraseasonal_wind < -threshold;   % Easterly wind events

% Connect nearby regions in time or space
wwe_regions = imclose(wwe_regions, se);
ewe_regions = imclose(ewe_regions, se);

% Label connected regions (4-connected components)
wwe_labels = bwlabel(wwe_regions, 4);
ewe_labels = bwlabel(ewe_regions, 4);

%% === Step 2: Define Seasonal and Longitudinal Bins ===
% Seasonal bins
% season_bins = {[11 12 1], [2 3 4], [5 6 7], [8 9 10]};
% season_names = {'NDJ', 'FMA', 'MJJ', 'ASO'};


season_bins = {1,2,3,4,5,6,7,8,9,10,11,12};
season_names = {'J', 'F','M','A', 'M','J','J', 'A','S','O','N','D'};

% Longitude bins exactly matching Figure 5
lon_edges = 120:20:280;             % Bin edges: 120, 140, ..., 280
lon_centers = 130:20:290;           % Tick labels: bin centers
n_bins = length(lon_centers);


% Initialize counters
season_wwe = zeros(1,12);  season_ewe = zeros(1,12);
location_wwe_bins = zeros(1,n_bins);  location_ewe_bins = zeros(1,n_bins);

%% === Step 3: Count Events by Season and Longitude Bin ===

% Westerly Wind Events (WWE)
for label = 1:max(wwe_labels(:))
    [x_idx,t_idx] = find(wwe_labels == label);
    if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
        mo = month(time_all(unique(t_idx)));  % Extract season
        for s = 1:12
            if any(ismember(mo, season_bins{s}))
                season_wwe(s) = season_wwe(s) + 1;
                break;
            end
        end
        lon_event = lon(unique(x_idx));
        event_center = mean(lon(x_idx));  % better than entire span
        for b = 1:n_bins
            if event_center >= lon_edges(b)-10 && event_center < lon_edges(b)+10
                location_wwe_bins(b) = location_wwe_bins(b) + 1;
                break;
            end
        end
    end
end


% Easterly Wind Events (EWE)
for label = 1:max(ewe_labels(:))
    [x_idx,t_idx] = find(ewe_labels == label);
    if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
        mo = month(time_all(unique(t_idx)));
        for s = 1:12
            if any(ismember(mo, season_bins{s}))
                season_ewe(s) = season_ewe(s) + 1;
                break;
            end
        end
        lon_event = lon(unique(x_idx));
        event_center = mean(lon(x_idx));
        for b = 1:n_bins
            if event_center >= lon_edges(b)-10 && event_center < lon_edges(b)+10
                location_ewe_bins(b) = location_ewe_bins(b) + 1;
                break;
            end
        end
    end
end

% Normalize to percentage
season_wwe_pct = 100 * season_wwe / sum(season_wwe);
season_ewe_pct = 100 * season_ewe / sum(season_ewe);
location_wwe_pct = 100 * location_wwe_bins / sum(location_wwe_bins);
location_ewe_pct = 100 * location_ewe_bins / sum(location_ewe_bins);

%% === Step 3: Plot Seasonality and Location Histograms ===
figure('position',[100 100 1000 600])

subplot(2,2,1)
bar(season_wwe_pct, 'FaceColor', [1 0.5 0])
set(gca, 'xticklabel', season_names, 'fontsize', 12)
ylabel('(%)'); title('(a) WWE Seasonality')

subplot(2,2,2)
bar(location_wwe_pct, 'FaceColor', [1 0.5 0])
set(gca, 'xtick', 1:n_bins, 'xticklabel', string(lon_edges), 'fontsize', 12)
xlabel('Longitude'); ylabel('(%)'); title('(b) WWE Location')

subplot(2,2,3)
bar(season_ewe_pct, 'FaceColor', [0.1 0.6 1])
set(gca, 'xticklabel', season_names, 'fontsize', 12)
ylabel('(%)'); title('(c) EWE Seasonality')

subplot(2,2,4)
bar(location_ewe_pct, 'FaceColor', [0.1 0.6 1])
set(gca, 'xtick', 1:n_bins, 'xticklabel', string(lon_edges), 'fontsize', 12)
xlabel('Longitude'); ylabel('(%)'); title('(d) EWE Location')




%% === Step 0: Select Subset Period for Display (e.g., 3 months)
% Pick time indices (e.g., 1 year)lon2
start_day = (1997-1979)*365+1;  % Choose your start date
end_day = (1997-1979)*365+1*365;    % Optional end date

% Convert time_all to datenum for indexing
time_nums = end_day-start_day+1;
%id_range = find(time_nums >= start_day & time_nums <= end_day);
id_range = start_day:end_day;

wind_slice = intraseasonal_wind(:, id_range);        % time × lon
wwe_mask = wwe_regions(:, id_range);                 % same size
ewe_mask = ewe_regions(:, id_range);

% Replace meshgrid line:
[XX, TT] = meshgrid(lon, (start_day:end_day)/365);

%% === Preparation ===
% use intraseasonal_wind (already selected Pacific region, left 120–right 280)
% transpose to [time × lon]
wind_field = intraseasonal_wind';  % [time × lon]
[ntime, nlon] = size(wind_field);
[LonGrid, TimeGrid] = meshgrid(lon, ttt);  % decimal year time vector

% extract event information
events = [];
for label = 1:max(wwe_labels(:))
    [x_idx,t_idx] = find(wwe_labels == label);
    if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
        t_evt = mean(ttt(unique(t_idx)));
        lon_evt = mean(lon(x_idx));
        strength = max(wind_field(t_idx, x_idx), [], 'all');
        events(end+1,:) = [t_evt, lon_evt, strength, 1];  % WWE
    end
end
for label = 1:max(ewe_labels(:))
    [x_idx,t_idx] = find(ewe_labels == label);
    if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
        t_evt = mean(ttt(unique(t_idx)));
        lon_evt = mean(lon(x_idx));
        strength = abs(min(wind_field(t_idx, x_idx), [], 'all'));
        events(end+1,:) = [t_evt, lon_evt, strength, -1];  % EWE
    end
end


%% === Step 0: Keep [lon × time] format ===
%u_band = bandpass(u_phy_obs, [1/100 1/20], 2);   % bandpass filtering
u_band = u_filt;
intraseasonal_wind = u_band(Left3_end:Right3_end,:);                    % already [lon × time]

lon_sub = lon;                                  % 122.5:2.5:280
ttt = 1979 + (0:size(u_band,2)-1)/365;
[~, ntime] = size(intraseasonal_wind);

%% === Step 1: Time windows ===
%years = [1995.5 1999.5; 2013.5 2017.5];
years = [2005.5 2008; 2013.5 2016.5];
titles = {'1997–1998', '2010'};


figure('position',[408.2,49.800000000000004,1047.2,1020.8])

for p = 1:2
    idx = ttt >= years(p,1) & ttt <= years(p,2);
    t_sub = ttt(idx);
    w_sub = intraseasonal_wind(:,idx);  % [lon × time]

    idx2 = idx(3*365+1:end);
    subplot(1,6,3*(p-1)+2)
    [TT, XX] = meshgrid(t_sub, lon_sub);
    contourf(XX, TT, w_sub, 20, 'LineStyle','none');  % background wind
    colormap('jet'); caxis([-10 10]); hold on;
    xlabel('Longitude'); ylabel('Year');
    colorbar('southoutside')
    title('(b) Intraseasonal Wind')
    
    set(gca,'fontsize',14)

    % === WWE: black contour lines for event regions ===
    for label = 1:max(wwe_labels(:))
        mask = (wwe_labels == label);
        [x_idx, t_idx] = find(mask);
        if isempty(x_idx), continue; end
        t_event = mean(ttt(t_idx));
        lon_event = mean(lon(x_idx));
        dur_t = length(unique(t_idx));
        span_x = range(lon(x_idx));
        if t_event < years(p,1) || t_event > years(p,2), continue; end
        if dur_t < min_duration || span_x < min_zonal_extent, continue; end

        mask_plot = double(mask(:, idx));
        if sum(mask_plot(:)) == 0, continue; end
        [Tgrid, Xgrid] = meshgrid(t_sub, lon_sub);
        scatter(lon_event, t_event, 5 + 5*max(w_sub(:)), 'MarkerFaceColor',[0.6 0 0.9], 'MarkerEdgeColor','k')
    end

    % === EWE: black dashed contour lines for event regions ===
    for label = 1:max(ewe_labels(:))
        mask = (ewe_labels == label);
        [x_idx, t_idx] = find(mask);
        if isempty(x_idx), continue; end
        t_event = mean(ttt(t_idx));
        lon_event = mean(lon(x_idx));
        dur_t = length(unique(t_idx));
        span_x = range(lon(x_idx));
        if t_event < years(p,1) || t_event > years(p,2), continue; end
        if dur_t < min_duration || span_x < min_zonal_extent, continue; end

        mask_plot = double(mask(:, idx));
        if sum(mask_plot(:)) == 0, continue; end
        [Tgrid, Xgrid] = meshgrid(t_sub, lon_sub);
        scatter(lon_event, t_event, 5 + 5*max(abs(w_sub(:))), 'MarkerFaceColor',[0.1 0.8 0.1], 'MarkerEdgeColor','k')
    end
    subplot(1,6,3*(p-1)+3)
    time2 = find(idx==1)-1095;
    [xx,yy] = meshgrid(X_sst(Left1_end:Right1_end),t_sub);
    contourf(xx,yy,sst_a_fine(time2,Left1_end:Right1_end),'LineStyle','none')%plot SST
    colormap('jet'); caxis([-3 3]); hold on;
    xlabel('Longitude'); ylabel('Year');
    set(gca,'fontsize',14)
    colorbar('southoutside')
    title('(c) SST')

    subplot(1,6,3*(p-1)+1)
    [TT, XX] = meshgrid(t_sub, lon_sub);
    contourf(XX,TT,MJO(Left3_end:Right3_end,idx),'LineStyle','none')%plot SST
    colormap('jet'); caxis([-0.7 0.7]); hold on;
    xlabel('Longitude'); ylabel('Year');
    set(gca,'fontsize',14)
    colorbar('southoutside')
    title('(a) MJO')
end

sgtitle('Intraseasonal Wind Events: 1997–1998 vs 2010','fontsize',18)