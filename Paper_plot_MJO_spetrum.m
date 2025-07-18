EigenSolver_rescale
st = 900*360;
lg = 1000;
Hbar = 22;

pos_x = [0.13, 0.57];
pos_y = [0.54,0.09];
wid = 0.26;
hig = 0.34;

% Extract time period data
u_phy_temp2 = u_phy(:,st+1:st+lg);
theta_phy_temp2 = theta_phy(:,st+1:st+lg);
K_phy_temp2 = K_phy(:,st+1:st+lg);
R_phy_temp2 = R_phy(:,st+1:st+lg);
Q_temp2 = Q(:,st+1:st+lg);
A_temp2 = A(:,st+1:st+lg);


jj = size(Q_temp2,2);
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

dim_t = 34;
xpoint = 5;
psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator

figure('Position',[744,353.8,772.2,596.2])

annotation('textbox',[0.408149749118661,0.897812814491781,0.18858818958819,0.04830593760483],'String','Observation','Fontsize',14,'LineStyle','none')
annotation('textbox',[0.396529396529397,0.449305602146931,0.206682206682207,0.048305937604831],'String','Model Simulation','Fontsize',14,'LineStyle','none')

for i = 1:2
    if i == 1
        % Use filtered K and R to calculate u
        variable = (K_phy_temp2-R_phy_temp2)*psi_0 + 1/sqrt(2)*R_phy_temp2*psi_2;
        variable = fftshift(fft(variable));
    elseif i == 2
        % A is already in Fourier space, use filtered A
        variable = A_temp2;
    end
    
    Spectrum_variable = zeros(size(variable));
    for j = 1:Na
        Spectrum_variable(j,:) = fftshift(fft(variable(j,:)'));
    end
    Spectrum_variable = sqrt(Spectrum_variable .* conj(Spectrum_variable))/length(variable(1,:));
    
    for j = 1:Na
        Spectrum_variable(j,:) = smooth(Spectrum_variable(j,:));
    end
    
    time_modes = round(length(variable(1,:))/10);
    xscale = -xpoint:xpoint;
    yscale = (1:time_modes)/(length(variable(1,:)));
    [xx_scale,yy_scale] = meshgrid(xscale,yscale);
    
    ax = subplot('Position',[pos_x(i),pos_y(2),wid,hig]);
    hold on
    contourf(ax,xx_scale, yy_scale, log(Spectrum_variable(zero_mode+xscale, zero_mode2 + (1:time_modes))'),40, 'linestyle','none');
    
    plot(0:5, Eig_Store(1:2,1:6)/2/pi/dim_t,'ko','linewidth',2)
    plot(-[0:5], -Eig_Store(3:4,1:6)/2/pi/dim_t,'ko','linewidth',2)
    xlim([-xpoint,xpoint])
    ylim([0,0.1])
    plot(ax,[-xpoint,xpoint],[1/30,1/30],'--k','linewidth',2)
    plot(ax,[-xpoint,xpoint],[1/90,1/90],'--k','linewidth',2)
    plot(ax,[0,0],[0,0.1],'--k','linewidth',2)
    box on
    set(gca,'FontSize',14)
    
    c = colorbar;
    set(c,'position',[pos_x(i)+0.28,pos_y(2),0.02,hig])
    set(ax,'Layer','top','Box','on','LineWidth',0.5)
end
colormap jet

load hgt_new_data
load uwnd_new_data
load Q_new_data
load olr_new_data

K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/sqrt(2);
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/sqrt(2) + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes);
A_3modes = olr_mode_0_rmseason_3modes;
Q_3modes = Q_mode_0_rmseason_3modes;

c = 0.05;
yy = -8:0.1:8;

% Atmosphere
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
 
u_phy2 = uwnd_mode_0_rmseason;
A_phy2 = olr_mode_0_rmseason;

u_phy2 = u_phy2(:,1:2000);
A_phy2 = A_phy2(:,1:2000);


jj = size(u_phy2,2);
if mod(jj,2) == 0
    zero_mode2 = (jj)/2+1;
else
    zero_mode2 = (jj+1)/2;
end

jj = size(u_phy2,1);
if mod(jj,2) == 0
    zero_mode = jj/2+1;
else
    zero_mode = (jj+1)/2;
end

% figure
for i = 1:2
    if i == 1
        % Use filtered observation data for u
        variable = u_phy2;
        variable = fftshift(fft(variable));
    elseif i == 2
        % Use filtered observation data for A
        variable = A_phy2;
        variable = fftshift(fft(variable));
    end
    Spectrum_variable = zeros(size(variable));
    for j = 1:size(u_phy2,1)    
        Spectrum_variable(j,:) = fftshift(fft(variable(j,:)'));
    end
    Spectrum_variable = sqrt(Spectrum_variable .* conj(Spectrum_variable))/length(variable(1,:));
    for j = 1:size(u_phy2,1)
        Spectrum_variable(j,:) = smooth(Spectrum_variable(j,:));
    end
    time_modes = round(length(variable(1,:))/10);
    xscale = -5:5;
    yscale = (1:time_modes)/(length(variable(1,:)));
    [xx_scale,yy_scale] = meshgrid(xscale,yscale);
    ax = subplot('Position',[pos_x(i),pos_y(1),wid,hig]);
    hold on
    contourf(ax,xx_scale, yy_scale, log(Spectrum_variable(zero_mode+xscale, zero_mode2 + (1:time_modes))'),40, 'linestyle','none');
    plot(ax,0:5, Eig_Store(1:2,1:6)/2/pi/dim_t,'ko','linewidth',2)
    plot(ax,-[0:5], -Eig_Store(3:4,1:6)/2/pi/dim_t,'ko','linewidth',2)
    xlim([-5,5])
    ylim([0,0.1])
    plot(ax,[-5,5],[1/30,1/30],'--k','linewidth',2)
    plot(ax,[-5,5],[1/90,1/90],'--k','linewidth',2)
    plot(ax,[0,0],[0,0.1],'--k','linewidth',2)
    box on
    set(gca,'FontSize',14)
    if i == 1
        title('(a) Zonal Velocity','Position',[0.00000955574066,0.117263781017675,0])
    elseif i == 2
        title('(b) Convectively Activity','Position',[0.000025703612552,0.117263781017675,0])
    end
    c = colorbar;
    set(c,'position',[pos_x(i)+0.28,pos_y(1),0.02,hig])
    set(ax,'Layer','top','Box','on','LineWidth',0.5)
end
colormap jet