steps_total = size(A_phy,2); % total time steps saved to disk
tt_decYr = (1:steps_total) / 360; % decimal-year vector
time_steps = steps_total;

% Fourier filtering parameters
time_high_period = 90; % high period cutoff (days) - upper limit of intraseasonal oscillation
time_low_period = 30; % low period cutoff (days) - lower limit of intraseasonal oscillation
max_wavenumber = 10; % maximum spatial wavenumber, keep wavenumbers from 0 to max_wavenumber (note: your original setting of 128 was too large)

% Extract data for the time range to be processed
u_slice = A_phy;

fprintf('Starting 2D Fourier spatiotemporal filtering...\n');
fprintf('Original data dimensions: [%d, %d]\n', size(u_slice,1), size(u_slice,2));

%% ---------- 1. 2D Fourier Transform -------------------------------------------
% Calculate 2D Fourier transform
fprintf('Calculating 2D Fourier transform...\n');

% Get actual data dimensions
[n_space, n_time] = size(u_slice);

% Calculate 2D FFT
u_fft2d = fft2(u_slice);

% Calculate frequency axes
fs_time = 1; % temporal sampling frequency (one point per day)
time_freqs = fs_time * (0:(n_time-1)) / n_time; % temporal frequency axis
spatial_wavenumbers = 0:(n_space-1); % spatial wavenumber axis

%% ---------- 2. Design 2D Filter -------------------------------------------
fprintf('Designing 2D filter...\n');

% Temporal frequency filter parameters
time_high_freq = 1/time_high_period; % low frequency corresponding to high period
time_low_freq = 1/time_low_period; % high frequency corresponding to low period

% Initialize filter
time_filter = zeros(n_space, n_time);
space_filter = zeros(n_space, n_time);

%% Temporal direction filter
for i = 1:n_space
    for j = 1:n_time
        f = time_freqs(j);

        % Handle frequency periodicity
        if (f >= time_high_freq && f <= time_low_freq) || ...
           (f >= (1-time_low_freq) && f <= (1-time_high_freq))
            time_filter(i, j) = 1;
        end
    end
end

%% Spatial direction filter - using smooth transition
space_transition_width = 3; % transition band width

for i = 1:n_space
    for j = 1:n_time
        k = spatial_wavenumbers(i);

        % Handle low wavenumber region
        if k <= max_wavenumber - space_transition_width
            space_filter(i, j) = 1; % complete pass
        elseif k <= max_wavenumber
            % smooth transition region
            weight = 0.5 * (1 + cos(pi * (k - max_wavenumber + space_transition_width) / space_transition_width));
            space_filter(i, j) = weight;

        % Handle high wavenumber region (considering FFT periodicity)
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

% Combine temporal and spatial filters
filter_2d = time_filter .* space_filter;

% Apply filter
A_fft2d_filtered = u_fft2d .* filter_2d;

%% ---------- 3. Inverse transform back to spatiotemporal domain -------------------------------------------
fprintf('Performing inverse Fourier transform...\n');

% Inverse transform
A_filtered = real(ifft2(A_fft2d_filtered));

% Verify dimensions
fprintf('Filtered data dimensions: [%d, %d]\n', size(A_filtered,1), size(A_filtered,2));
fprintf('Original data dimensions: [%d, %d]\n', size(A_phy,1), size(A_phy,2));

% Ensure output dimensions match input
if ~isequal(size(A_filtered), size(A_phy))
    error('Filtered dimensions do not match! Please check code');
end

fprintf('2D filtering completed!\n');

Starting_year = 495;
Ending_year = 525;

LL = 10;

dim_t = 34;
dt = 8/24/dim_t/10; % time step, every dt*dim_t hours
gap = 30;

LL2 = floor((Starting_year*360+1)/dim_t/(dt*gap)); % Starting time of plotting
LL3 = floor(Ending_year*360/dim_t/(dt*gap)); % Final time of plotting
ll = length(dt*LL2:dt: dt*LL3);
idx = LL2:LL:LL2+ll-1;
idx2 = idx;



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

T4_store_a  = movmean(T4_store_a, 30);
T3_store_a  = movmean(T3_store_a, 30);
T34_store_a = movmean(T34_store_a, [1 90]);

temp_llon = 0:360/(Na):360-360/(Na);

temp_x = (-Na+No-1:No) * dx; % grid points in the x axis
temp_x2 = temp_llon;

[temp_xx,temp_yy] = meshgrid(temp_x*dim_x/1000,idx2/360);
[temp_xx2,temp_yy2] = meshgrid(temp_x2,idx2/360);



% Get the length of the time dimension
lt = size(A_phy, 2);

% Calculate the zero mode (center of the frequency domain)
zero_mode = floor(lt/2) + 1;

temp_A_phy = A_filtered([No+1:Na,1:No],:);

olr_mode_0_rmseason = temp_A_phy;

% Parameters
[~, ntime] = size(olr_mode_0_rmseason);
time = 1:ntime;

olr_filtered = olr_mode_0_rmseason;


olr_amplitude = sqrt(movmean(olr_filtered.^2, 360,2))/10;

[xx,yy] = meshgrid(temp_llon,idx/360);


figure
t = tiledlayout(1,2);
ax1 = axes(t);
plot(ax1,mean(olr_amplitude(72+(1:No/2),idx),1), idx/360,'-r')
ax2 = axes(t);
plot(ax2,T34_store_a(idx), idx/360,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax1.YLim = [idx(1)/360,idx(end)/360];
ax2.YLim = [idx(1)/360,idx(end)/360];

figure
colormap jet
% First subplot remains the same
subplot(1,4,3);
contourf(xx, yy, olr_amplitude(:,idx)','LineStyle','none');
title('Interannual Amplitude of Intraseasonal OLR');
set(gca, 'XTick', 20:60:320);
set(gca, 'XTickLabel', {'60°W', '0°E', '60°E', '120°E', '180°', '120°W'});
xlabel('Longitude');
ylabel('Time (years)');
clim([0,0.15])
colorbar;

load Obs_fine_2020.mat
load olr_new_data

time2 = 1:10950;
time3 = 1095+time2;

Left1_end = (120+2.5) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST
Middle1_end = 200 / 2.5 * 10; % Middle of Pacific ocean 160W for SST
Left1_C = 160/2.5*10;
Right1_C = 210/2.5*10;
Left1_C2 = 190/2.5*10;
Right1_C2 = 240/2.5*10;

Left3_end = (120+2.5) / 2.5; % Pacific ocean left boundary 120E for wind bursts
Right3_end = 280 / 2.5; % Pacific ocean right boundary 80W for wind bursts
Middle3_end = 200 / 2.5; % Middle of Pacific ocean 160W for wind bursts
Left3_1 = (30+2.5)/2.5;

T34_new = mean(sst_a_fine(time2,Left1_C2:Right1_C2),2);
T4_new = mean(sst_a_fine(time2,Left1_C:Right1_C),2);
T3_new = mean(sst_a_fine(time2,Middle1_end:Right1_end),2);
T34_new = movmean(T34_new, [1,90]);
%T34_new = T34_new(1:30:end);

olr_mode_0_rmseason = olr_mode_0_rmseason(:,time3);
olr_mode_0_rmseason_3modes = olr_mode_0_rmseason_3modes(:,time3);
% Parameters
[~, ntime] = size(olr_mode_0_rmseason);
time = 1:ntime;


bd1 = -90; % Lower bound (100-day period)
bd2 = -30;  % Upper bound (20-day period)

% Get the length of the time dimension
lt = size(olr_mode_0_rmseason, 2);

% Calculate the zero mode (center of the frequency domain)
zero_mode = floor(lt/2) + 1;


olr_filtered = zeros(size(olr_mode_0_rmseason));
for i = 1:144
    Ha_temporal = fftshift(fft(olr_mode_0_rmseason(i,:)));
    Ha_temporal2 = Ha_temporal * 0;
    Ha_temporal2(zero_mode+round(lt/bd2):zero_mode+round(lt/bd1)) = Ha_temporal(zero_mode+round(lt/bd2):zero_mode+round(lt/bd1));
    Ha_temporal2(zero_mode-round(lt/bd1):zero_mode-round(lt/bd2)) = Ha_temporal(zero_mode-round(lt/bd1):zero_mode-round(lt/bd2));
    olr_filtered(i,:) = 2*real(ifft(ifftshift(Ha_temporal2)));
end


olr_filtered = olr_mode_0_rmseason_3modes;

% 2. Calculate interannual amplitude
olr_amplitude = sqrt(movmean(olr_filtered.^2, 365,2));

% 4. Calculate running correlation between PC1 and PC3
window = 400;

% Create the second subplot with dual x-axes
temp_lon = [(282.5+2.5)/2.5:(357.5+2.5)/2.5,1:(280+2.5)/2.5];
[xx,yy] = meshgrid(lon,Y(time2));

colormap jet
% First subplot remains the same
subplot(1,4,1);
contourf(xx, yy, olr_amplitude(temp_lon,:)','LineStyle','none');
title('Interannual Amplitude of Intraseasonal OLR');
set(gca, 'XTick', 20:60:320);
set(gca, 'XTickLabel', {'60°W', '0°E', '60°E', '120°E', '180°', '120°W'});
xlabel('Longitude');
ylabel('Time (years)');
colorbar;
clim([0.05,0.15])

figure
t = tiledlayout(1,2);
ax1 = axes(t);
plot(ax1,mean(olr_amplitude(Left3_end:Middle3_end,:),1), Y(time2),'-r')
ax2 = axes(t);
plot(ax2,T34_new, Y(time2),'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax1.YLim = [Y(time2(1)),Y(time2(end))];
ax2.YLim = [Y(time2(1)),Y(time2(end))];