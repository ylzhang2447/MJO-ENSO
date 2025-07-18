min_duration_day = 3; % minimum duration (days)
min_span_deg = 10; % minimum zonal span (degrees)
dim_t = 34; % in model: 8-h is divided into dim_t integration steps
gap = 30; % write to file every gap dt steps

% data dimensions
Na = 128; % number of longitude points (spatial direction)
No = 56; % number of latitude points (if applicable)
total_steps = 360*2000; % total time steps

% Fourier filtering parameters
time_high_period = 90; % high period cutoff (days) - upper limit of intraseasonal oscillation
time_low_period = 30; % low period cutoff (days) - lower limit of intraseasonal oscillation
max_wavenumber = 10; % maximum spatial wavenumber, keep wavenumbers from 0 to max_wavenumber

%% ---------- 0. Derived grid and time axis -----------------------------------------
%lon = linspace(30, 280, Na); % align with zonal grid (Na points)
lon = 120+1/No:160/No:280;

% create time axis (in days)
days_per_year = 360; % based on 360-day calendar
total_years = total_steps / days_per_year;
time_axis = linspace(1, total_years, total_steps);

% determine time range to process
selected_time = time_axis;
time_steps2 = length(time_axis);
time_steps = total_steps;

% extract data for the time range to be processed
u_slice = 1.8*u_phy;

% get actual data dimensions
[n_space, n_time] = size(u_slice);

% calculate 2D FFT
u_fft2d = fft2(u_slice);

% calculate frequency axes
fs_time = 1; % temporal sampling frequency (one point per day)
time_freqs = fs_time * (0:(n_time-1)) / n_time; % temporal frequency axis
spatial_wavenumbers = 0:(n_space-1); % spatial wavenumber axis

time_high_freq = 1/time_high_period; 
time_low_freq = 1/time_low_period;

time_filter = zeros(n_space, n_time);
space_filter = zeros(n_space, n_time);

for i = 1:n_space
    for j = 1:n_time
        f = time_freqs(j);
        if (f >= time_high_freq && f <= time_low_freq) || ...
           (f >= (1-time_low_freq) && f <= (1-time_high_freq))
            time_filter(i, j) = 1;
        end
    end
end

for i = 1:n_space
    k = spatial_wavenumbers(i);
    
    % use sharp cutoff
    if k <= max_wavenumber || k >= (n_space - max_wavenumber)
        space_filter(i, :) = 1;
    else
        space_filter(i, :) = 0;
    end
end

% combine temporal and spatial filters
filter_2d = time_filter .* space_filter;

% apply filter
u_fft2d_filtered = u_fft2d .* filter_2d;

% inverse transform
u_filtered = real(ifft2(u_fft2d_filtered));

u_filt = u_filtered;



%% === Step 0: Prepare Intraseasonal Wind Data and Time ===
%u_filt = bandpass(u_phy_obs, [1/100 1/20], 2);  % filter along time dim (dim=2)
%u_filt = u_phy;
time_total = 2000;
%u_filt = movmean(u_phy_obs, 30, 2, 'omitnan', 'Endpoints', 'fill');
intraseasonal_wind = u_filt(1:No,1:time_total*360);  % [time × lon]
ttt = [1:360*time_total]/360;              % decimal year time vector
[nlon,ntime] = size(intraseasonal_wind);

%u_std = std(u_phy_obs(:), 'omitnan');
threshold = 5;  % e.g., 95th percentile filtering
% Event detection parameters
%threshold = 2;            % Wind stress threshold (N/m²)
min_duration = 3;         % Minimum duration (days)
min_zonal_extent = 10;     % Minimum zonal span (degrees)
se = strel('rectangle', [1, 3]);  % Structuring element for morphological closing

% Convert decimal year to datetime
%time_all = datetime(floor(ttt), 1, 1) + calmonths(round(12 * (ttt - floor(ttt))));


%% === Step 1: Binary Event Masking and Morphological Closing ===
wwe_regions = intraseasonal_wind > threshold;    % Westerly wind events
ewe_regions = intraseasonal_wind < -threshold;   % Easterly wind events

% Connect nearby regions in time or space
wwe_regions = imclose(wwe_regions, se);
ewe_regions = imclose(ewe_regions, se);

% Label connected regions (4-connected components)
wwe_labels = bwlabel(wwe_regions, 4);
ewe_labels = bwlabel(ewe_regions, 4);


n_segments = 45;
segment_years = 44;
n_seasons = 12;
% season_bins = {[11 12 1], [2 3 4], [5 6 7], [8 9 10]};
% season_names = {'NDJ', 'FMA', 'MJJ', 'ASO'};
season_bins = {1,2,3,4,5,6,7,8,9,10,11,12};
season_names = {'J', 'F','M','A', 'M','J','J', 'A','S','O','N','D'};

lon_edges = 120:20:280;
lon_centers = 130:20:290;
n_bins = length(lon_centers);

% Store per-segment stats
season_wwe_all = zeros(n_segments, n_seasons);
season_ewe_all = zeros(n_segments, n_seasons);
location_wwe_all = zeros(n_segments, n_bins);
location_ewe_all = zeros(n_segments, n_bins);

% ==================== Loop over segments =====================
for s = 1:n_segments
    idx_start = (s-1)*360*segment_years + 1;
    idx_end   = s*360*segment_years;

    wwe_seg = wwe_labels(:, idx_start:idx_end);
    ewe_seg = ewe_labels(:, idx_start:idx_end);

    for label = 1:max(wwe_seg(:))
        [x_idx, t_idx] = find(wwe_seg == label);
        if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
            mo = mod(floor((unique(t_idx)-1)/30), 12) + 1;
            for k = 1:n_seasons
                if any(ismember(mo, season_bins{k}))
                    season_wwe_all(s,k) = season_wwe_all(s,k) + 1;
                    break;
                end
            end
            event_center = mean(lon(x_idx));
            for b = 1:n_bins
                if event_center >= lon_edges(b)-10 && event_center < lon_edges(b)+10
                    location_wwe_all(s,b) = location_wwe_all(s,b) + 1;
                    break;
                end
            end
        end
    end

    for label = 1:max(ewe_seg(:))
        [x_idx, t_idx] = find(ewe_seg == label);
        if length(unique(t_idx)) >= min_duration && range(lon(x_idx)) >= min_zonal_extent
            mo = mod(floor((unique(t_idx)-1)/30), 12) + 1;
            for k = 1:n_seasons
                if any(ismember(mo, season_bins{k}))
                    season_ewe_all(s,k) = season_ewe_all(s,k) + 1;
                    break;
                end
            end
            event_center = mean(lon(x_idx));
            for b = 1:n_bins
                if event_center >= lon_edges(b)-10 && event_center < lon_edges(b)+10
                    location_ewe_all(s,b) = location_ewe_all(s,b) + 1;
                    break;
                end
            end
        end
    end
end

% ==================== Normalize to percentage =====================
season_wwe_pct_all = 100 * season_wwe_all ./ sum(season_wwe_all,2);
season_ewe_pct_all = 100 * season_ewe_all ./ sum(season_ewe_all,2);
location_wwe_pct_all = 100 * location_wwe_all ./ sum(location_wwe_all,2);
location_ewe_pct_all = 100 * location_ewe_all ./ sum(location_ewe_all,2);

% ==================== Compute mean and std =====================
mean_season_wwe = mean(season_wwe_pct_all,1);
std_season_wwe = std(season_wwe_pct_all,0,1);
mean_season_ewe = mean(season_ewe_pct_all,1);
std_season_ewe = std(season_ewe_pct_all,0,1);

mean_loc_wwe = mean(location_wwe_pct_all,1);
std_loc_wwe = std(location_wwe_pct_all,0,1);
mean_loc_ewe = mean(location_ewe_pct_all,1);
std_loc_ewe = std(location_ewe_pct_all,0,1);

% ==================== Plotting =====================
figure('Position',[100 100 1200 600])

subplot(2,2,1)
bar(mean_season_wwe, 'FaceColor', [1 0.5 0]); hold on
errorbar(1:n_seasons, mean_season_wwe, std_season_wwe, 'k.', 'LineWidth', 1.2)
set(gca,'xticklabel',season_names, 'fontsize',12)
ylabel('WWE (%)'); title('(a) WWE Seasonality')

subplot(2,2,2)
bar(mean_loc_wwe, 'FaceColor', [1 0.5 0]); hold on
errorbar(1:n_bins, mean_loc_wwe, std_loc_wwe, 'k.', 'LineWidth', 1.2)
set(gca,'xtick',1:n_bins,'xticklabel',string(lon_centers-10),'fontsize',12)
xlabel('Longitude'); ylabel('WWE (%)'); title('(b) WWE Location')

subplot(2,2,3)
bar(mean_season_ewe, 'FaceColor', [0.1 0.6 1]); hold on
errorbar(1:n_seasons, mean_season_ewe, std_season_ewe, 'k.', 'LineWidth', 1.2)
set(gca,'xticklabel',season_names, 'fontsize',12)
ylabel('EWE (%)'); title('(c) EWE Seasonality')

subplot(2,2,4)
bar(mean_loc_ewe, 'FaceColor', [0.1 0.6 1]); hold on
errorbar(1:n_bins, mean_loc_ewe, std_loc_ewe, 'k.', 'LineWidth', 1.2)
set(gca,'xtick',1:n_bins,'xticklabel',string(lon_centers-10),'fontsize',12)
xlabel('Longitude'); ylabel('EWE (%)'); title('(d) EWE Location')
