%% Time-filtering version of u_phy Fourier filtering code (30-90 days)

steps_total = size(u_phy,2); % Total integration steps
tt_decYr = (1:steps_total) / 360; % Decimal-year vector

%% ======================= User-defined parameters ========================
yr0 = 155;
yr1 = 158;
threshold = 5; % m s-1, wind speed threshold (±threshold)
min_duration_day = 3; % Minimum duration (days)
min_span_deg = 10; % Minimum meridional span (degrees)
dim_t = 34; % In model: 8-h is divided into dim_t integration steps
gap = 30; % Write file every gap dt steps

% Data dimensions
Na = 128; % Number of longitude points (spatial direction)
No = 56; % Number of latitude points (if applicable)

% Time filtering parameters
time_high_period = 90; % High period cutoff (days) - upper limit of intraseasonal oscillation
time_low_period = 30; % Low period cutoff (days) - lower limit of intraseasonal oscillation

%% ---------- 0. Derived grid and time axis ----------------------------------
lon = linspace(30, 280, Na); % Aligned with meridional grid (Na points)
dt = 1/(gap*dim_t); % Single integration step (days)

fprintf('Original data dimensions: [%d, %d]\n', size(u_phy));

% Get actual data dimensions
[n_space, n_time] = size(u_phy);

fprintf('Starting time filtering (retaining %d-%d day periods)...\n', time_low_period, time_high_period);

%% ---------- 1. Time filtering (along time direction) ---------------------

% Calculate time zero mode position
if mod(n_time, 2) == 0
    zero_mode_time = n_time/2 + 1;
else
    zero_mode_time = (n_time + 1)/2;
end

% Calculate frequency boundaries (similar to reference code bd1, bd2)
% Note: negative sign used here due to reference code convention
bd1 = time_high_period;  % Corresponds to 90 days
bd2 = time_low_period;   % Corresponds to 30 days

fprintf('Time filtering boundaries: bd1=%.2f, bd2=%.2f\n', bd1, bd2);
fprintf('Time zero mode position: %d\n', zero_mode_time);

% Initialize time-filtered data
u_filtered = zeros(n_space, n_time);

% Apply time filtering to each spatial point
for i = 1:n_space
    % FFT of current spatial point's time series
    u_temporal = fftshift(fft(u_phy(i, :)));
    
    % Create time filter
    u_temporal_filtered = u_temporal * 0;  % Initialize to zero
    
    % Retain components in specified frequency range
    freq_start = zero_mode_time + round(n_time/bd1);
    freq_end = zero_mode_time + round(n_time/bd2);
    
    % Ensure indices are within valid range
    freq_start = max(1, min(freq_start, n_time));
    freq_end = max(1, min(freq_end, n_time));
    
    if freq_start <= freq_end
        u_temporal_filtered(freq_start:freq_end) = u_temporal(freq_start:freq_end);
    end
    
    % For negative frequency part (symmetry)
    neg_freq_start = zero_mode_time - round(n_time/bd2);
    neg_freq_end = zero_mode_time - round(n_time/bd1);
    
    neg_freq_start = max(1, min(neg_freq_start, n_time));
    neg_freq_end = max(1, min(neg_freq_end, n_time));
    
    if neg_freq_start <= neg_freq_end
        u_temporal_filtered(neg_freq_start:neg_freq_end) = u_temporal(neg_freq_start:neg_freq_end);
    end
    
    % Inverse transform back to time domain
    u_filtered(i, :) = real(ifft(ifftshift(u_temporal_filtered)));
end

fprintf('Time filtering completed!\n');

% Verify dimensions
fprintf('Filtered data dimensions: [%d, %d]\n', size(u_filtered,1), size(u_filtered,2));
fprintf('Original data dimensions: [%d, %d]\n', size(u_phy,1), size(u_phy,2));

% Use filtered data
u_new = 1.8*u_filtered;

%% ---------- 3. Subsequent analysis code remains unchanged ----------------
% Extract specified time period
mask_time  = (tt_decYr >= yr0) & (tt_decYr <= yr1);
u_sub      = u_new(1:No, mask_time);      % [lon × time_sub]
T_sub      = Hov_T_a(1:No, mask_time);    % [lon × time_sub]
tt_sub     = tt_decYr(mask_time);         % Synchronized time axis
MJO_sub      = MJO(1:No, mask_time);      % [lon × time_sub]

% Threshold masking + morphological closing operation
se        = strel('rectangle',[1 3]);     
wwe_bin   = imclose(u_sub >  threshold, se);
ewe_bin   = imclose(u_sub < -threshold, se);

% Connected component filtering
dt_day = 1;                               
wwe_lbl = filterByExtent(bwlabel(wwe_bin,4), lon, dt_day, ...
                         min_duration_day, min_span_deg);
ewe_lbl = filterByExtent(bwlabel(ewe_bin,4), lon, dt_day, ...
                         min_duration_day, min_span_deg);

% Plotting section remains unchanged...
temp_lon = (-32+1:No) * 360/Na + 120;     
sub_lon  = temp_lon(end-No+1:end);        
[XX,YY]  = meshgrid(temp_lon, tt_sub);    

figure
subplot(1,3,2)
contourf(XX(:,end-No+1:end), YY(:,end-No+1:end), u_sub', ...
         20, 'LineStyle','none');
colormap(jet);  caxis([-10 10]); colorbar('southoutside')
xlabel('Longitude');  ylabel('Year')
title('Intraseasonal Wind (30-90 day filtered)')
set(gca,'fontsize',14); hold on

% Overlay WWE (Westerly Wind Events)
for lab = 1:max(wwe_lbl(:))
    mask = (wwe_lbl==lab);  if ~any(mask(:)), continue, end
    [x,t]  = find(mask);                
    scatter(mean(sub_lon(x)), mean(tt_sub(t)), 5 + 5*max(abs(u_sub(:))), ...
            'MarkerFaceColor',[0.6 0 0.9],'MarkerEdgeColor','k');
end

% Overlay EWE (Easterly Wind Events)
for lab = 1:max(ewe_lbl(:))
    mask = (ewe_lbl==lab);  if ~any(mask(:)), continue, end
    [x,t]  = find(mask);
    scatter(mean(sub_lon(x)), mean(tt_sub(t)), 5 + 5*max(abs(u_sub(:))), ...
            'MarkerFaceColor',[0.1 0.8 0.1],'MarkerEdgeColor','k');
end

subplot(1,3,3)
contourf(XX(:,end-No+1:end), YY(:,end-No+1:end), T_sub', ...
         20, 'LineStyle','none');
colormap(jet);  caxis([-3.5 3.5]); colorbar('southoutside')
xlabel('Longitude');  ylabel('Year')
title('SST')
set(gca,'fontsize',14); hold on

subplot(1,3,1)
contourf(XX(:,end-No+1:end), YY(:,end-No+1:end), MJO_sub'/2.5, ...
         20, 'LineStyle','none');
colormap(jet);  caxis([-0.7 0.7]); 
colorbar('southoutside')
xlabel('Longitude');  ylabel('Year')
title('MJO')
set(gca,'fontsize',14); hold on

% Also overlay WWE and EWE on MJO plot
for lab = 1:max(wwe_lbl(:))
    mask = (wwe_lbl==lab);  if ~any(mask(:)), continue, end
    [x,t]  = find(mask);                
    scatter(mean(sub_lon(x)), mean(tt_sub(t)), 5 + 5*max(abs(u_sub(:))), ...
            'MarkerFaceColor',[0.6 0 0.9],'MarkerEdgeColor','k');
end

for lab = 1:max(ewe_lbl(:))
    mask = (ewe_lbl==lab);  if ~any(mask(:)), continue, end
    [x,t]  = find(mask);
    scatter(mean(sub_lon(x)), mean(tt_sub(t)), 5 + 5*max(abs(u_sub(:))), ...
            'MarkerFaceColor',[0.1 0.8 0.1],'MarkerEdgeColor','k');
end

%% ---------- 4. Connected component filtering function definition ----------
function lbl_out = filterByExtent(lbl_in, lon_vec, dt_day, ...
                                  min_dur_day, min_span_deg)
    if all(lbl_in(:)==0), lbl_out = lbl_in; return, end
    deg_per_grid = mean(diff(lon_vec));               
    min_span_pts = ceil(min_span_deg / deg_per_grid); 
    min_dur_pts  = ceil(min_dur_day  / dt_day);       

    stats = regionprops(lbl_in, 'PixelIdxList','BoundingBox');
    for k = 1:numel(stats)
        bb   = stats(k).BoundingBox;     
        span = bb(3);  dur = bb(4);
        if span < min_span_pts || dur < min_dur_pts
            lbl_in(stats(k).PixelIdxList) = 0;   
        end
    end
    lbl_out = bwlabel(lbl_in>0,4);       
end