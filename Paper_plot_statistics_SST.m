%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison of the model statistics with observations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the observational data nino 3, 4 and 3.4
data=importdata('.\nino3.txt');
nino3=reshape(data(:,2:end)',[],1);
nino3=nino3((1980-1870)*12+1:end-12);
%nino3=nino3((1951-1870)*12+1:end-12);
nino3_seasonal = reshape(nino3,12,[]);
nino3_seasonal = var(nino3_seasonal');

data=importdata('.\nino4.txt');
nino4=reshape(data(:,2:end)',[],1);
nino4=nino4((1980-1870)*12+1:end-12);
nino4_seasonal = reshape(nino4,12,[]);
nino4_seasonal = var(nino4_seasonal');

data=importdata('.\nino34.txt');
nino34=reshape(data(:,2:end)',[],1);
nino34=nino34((1980-1870)*12+1:end-12);
nino34_seasonal = reshape(nino34,12,[]);
nino34_seasonal = var(nino34_seasonal');

figure('Position',[494.6,350.6,1464,599.2])
% Show PDFs and the seasonal variation of the variance in different Nino
% regions, comparison with observations
posx = [0.04,0.2,0.36,0.54,0.7,0.86];
posy = [0.6,0.1];
posw = 0.13; posh = 0.3;

% nino 4
subplot('Position',[posx(1),posy(1),posw,posh])
%subplot(2,3,1)
[fi,xx]=ksdensity(T4_store_a); fi = smooth(fi);
[fi_obs,xx] = ksdensity(nino4,xx);
L = length(nino4);
times = floor(length(T4_store_a)/L);
PDF_T_C_model=zeros(times,length(fi_obs));
for i=1:times
    [PDF_T_C_model(i,:), xx] = ksdensity(T4_store_a(1+(i-1)*L:L+(i-1)*L)-mean(T4_store_a),xx);
end
PDF_T_C_model_m=nanmean(PDF_T_C_model,1);
PDF_T_C_model_std=nanstd(PDF_T_C_model,0,1);
facealpha = 0.6;
hold on 
hh1 = fill([xx fliplr(xx)],...
    [PDF_T_C_model_m-PDF_T_C_model_std fliplr(PDF_T_C_model_m+PDF_T_C_model_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(xx,fi_obs,'r','linewidth',2)
plot(xx,PDF_T_C_model_m,'b','linewidth',2)
set(gca,'fontsize',14)
title('Nino 4')
xlim([-3,4.2])
ylim([0,1.1])
box on
%legend('Model','Obs')
ylabel('probability')
ylim([0,0.9])
xlabel('^oC')
% nino 3.4
%subplot(2,3,2)
subplot('Position',[posx(2),posy(1),posw,posh])
[fi,xx]=ksdensity(T34_store_a);fi = smooth(fi);
[fi_obs,xx] = ksdensity(nino34,xx);
L = length(nino4);
times = floor(length(T34_store_a)/L);
PDF_T_34_model=zeros(times,length(fi_obs));
for i=1:times
    [PDF_T_34_model(i,:), xx] = ksdensity(T34_store_a(1+(i-1)*L:L+(i-1)*L)-mean(T34_store_a),xx);
end
PDF_T_34_model_m = nanmean(PDF_T_34_model,1);
PDF_T_34_model_std = nanstd(PDF_T_34_model,0,1);
facealpha = 0.6;
hold on
hh1 = fill([xx fliplr(xx)],...
    [PDF_T_34_model_m-PDF_T_34_model_std fliplr(PDF_T_34_model_m+PDF_T_34_model_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(xx,fi_obs,'r','linewidth',2)
plot(xx,PDF_T_34_model_m,'b','linewidth',2)
set(gca,'fontsize',14)
title({['(a) PDFs'], ['Nino 3.4']})
xlim([-3,4.2])
ylim([0,1.1])
box on
ylim([0,0.9])
xlabel('^oC')
% nino 3
subplot('Position',[posx(3),posy(1),posw,posh])
%subplot(2,3,3)
[fi,xx]=ksdensity(T3_store_a);fi = smooth(fi);
[fi_obs,xx] = ksdensity(nino3,xx);
hold on
PDF_T_E_model=zeros(times,length(fi_obs));
for i=1:times
    [PDF_T_E_model(i,:), xx] = ksdensity(T3_store_a(1+(i-1)*L:L+(i-1)*L)-mean(T3_store_a),xx);
end
PDF_T_E_model_m = nanmean(PDF_T_E_model,1);
PDF_T_E_model_std = nanstd(PDF_T_E_model,0,1);
facealpha = 0.6;
hold on
hh1 = fill([xx fliplr(xx)],...
    [PDF_T_E_model_m-PDF_T_E_model_std fliplr(PDF_T_E_model_m+PDF_T_E_model_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(xx,fi_obs,'r','linewidth',2)
plot(xx,PDF_T_E_model_m,'b','linewidth',2)
set(gca,'fontsize',14)
title('Nino 3')
xlim([-3,4.2])
ylim([0,1.1])
box on
ylim([0,0.9])
xlabel('^oC')

% plot seasonal variation of the variance
T_3_seasonal=zeros(times,12);
T_4_seasonal=zeros(times,12);
T_34_seasonal=zeros(times,12);
for i=1:times
    temp=reshape(T3_store_a(1+(i-1)*L:L+(i-1)*L),12,[]);
    T_3_seasonal(i,:)=var(temp');
    temp=reshape(T4_store_a(1+(i-1)*L:L+(i-1)*L),12,[]);
    T_4_seasonal(i,:)=var(temp');
    temp=reshape(T34_store_a(1+(i-1)*L:L+(i-1)*L),12,[]);
    T_34_seasonal(i,:)=var(temp');
end
T_3_seasonal_m=mean(T_3_seasonal,1);
T_3_seasonal_std=std(T_3_seasonal,0,1);
T_4_seasonal_m=mean(T_4_seasonal,1);
T_4_seasonal_std=std(T_4_seasonal,0,1);
T_34_seasonal_m=mean(T_34_seasonal,1);
T_34_seasonal_std=std(T_34_seasonal,0,1);

subplot('Position',[posx(3),posy(2),posw,posh])
%subplot(2,3,6)
hold on 
hh1 = fill([1:12 fliplr(1:12)],...
    [T_3_seasonal_m-T_3_seasonal_std fliplr(T_3_seasonal_m+T_3_seasonal_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(1:12,T_3_seasonal_m,'b','linewidth',2)
plot(1:12,nino3_seasonal,'r','linewidth',2)
set(gca,'fontsize',14)
set(gca,'xlim',[1 12]);
title('Nino 3')
xlabel('Calendar month')
box on
%ylabel('^oC^2')
ylim([0,2.2])

subplot('Position',[posx(2),posy(2),posw,posh])
%subplot(2,3,5)
hold on
hh1 = fill([1:12 fliplr(1:12)],...
    [T_34_seasonal_m-T_34_seasonal_std fliplr(T_34_seasonal_m+T_34_seasonal_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(1:12,T_34_seasonal_m,'b','linewidth',2)
plot(1:12,nino34_seasonal,'r','linewidth',2)
set(gca,'fontsize',14)
set(gca,'xlim',[1 12]);
title({['(c) Seasonal variances'], ['Nino 3.4']})
xlabel('Calendar month')
box on
ylim([0,2.2])

subplot('Position',[posx(1),posy(2),posw,posh])
%subplot(2,3,4)
hold on 
hh1 = fill([1:12 fliplr(1:12)],...
    [T_4_seasonal_m-T_4_seasonal_std fliplr(T_4_seasonal_m+T_4_seasonal_std)],'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot(1:12,T_4_seasonal_m,'b','linewidth',2)
plot(1:12,nino4_seasonal,'r','linewidth',2)
set(gca,'fontsize',14)
set(gca,'xlim',[1 12]);
title('Nino 4')
xlabel('Calendar month')
box on
ylim([0,2.2])
ylabel('(^oC)^2')

%figure
% plot the spectrum
num_piece = floor(length(T3_store_a)/480)-1;
save_spec = zeros(257, num_piece);
for i = 1:3
    if i == 1
        index1all = T4_store_a;
        index2 = nino4;
    elseif i == 3
        index1all = T34_store_a;
        index2 = nino34;
    elseif i == 2
        index1all = T3_store_a;
        index2 = nino3;
    end
    subplot('Position',[posx(3+i),posy(1),posw,posh])
    %subplot(2,3,i)
    for j = 1:num_piece
        ts = 1*12:1*12:480*12;
        Fs = 1*12;
        index1=index1all((j)*480+[1:480]);
        L = length(ts);
        NFFT = 2^nextpow2(L); 
        Y_y1 = fft(index1,NFFT)/L;
        f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
        tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
        save_spec(:,j) = tpp_y;
    end
    mean_model= mean(save_spec');
    std_model = std(save_spec');
    hold on
    plot(1./f_y1(end:-1:1) ,mean_model(end:-1:1),'b','linewidth',2) 
    low = mean_model-std_model;low(low<1e-10)=1e-10;
    x_axis = [f_y1(end:-1:1), fliplr(f_y1(end:-1:1))]; x_axis(x_axis==0)=0.00001;x_axis = 1./x_axis; 
    y_axis = [low(end:-1:1) fliplr(mean_model(end:-1:1)+std_model(end:-1:1))];  
    facealpha=0.6;
    hh1 = fill(x_axis, y_axis,'b','facealpha',facealpha);
    hh1.EdgeColor = 'none';
    set(gca,'xscale','log')
    xlim([0.8,10]);
    set(gca,'fontsize',14)
    box on
    xlabel('Year')
    set(gca,'xTick',[1:6,8,10])
    set(gca,'xTicklabel',  [1:6,8,10]);
    hold on
    ts = 1*12:1*12:480*12; 
    Fs = 1*12; 
    L = length(ts);
    NFFT = 2^nextpow2(L); 
    Y_y1 = fft(index2,NFFT)/L;
    f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
    tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
    hh2=plot(1./f_y1(end:-1:1),tpp_y(end:-1:1),'r','linewidth',2);
    temp_y = 2*abs(Y_y1(1:NFFT/2+1));
    tp_y = temp_y; 
    if i == 1
        title('Nino 4')
        ylabel('power')
        legend([hh1,hh2],'Model','Obs')
    elseif i == 2
        title({['(b) Spectrums'], ['Nino 3.4']})
    else
        title('Nino 3')
    end
    ylim([0,0.8])
    % ylim([1e-3,1])

end

% plot the autocorrelation function (ACF)

load('obs_data');
ACF_T_E_model = autocorr(T3_store_a,60);
ACF_T_C_model = autocorr(T4_store_a,60);
ACF_T_CE_model = autocorr(T34_store_a,60);

ACF_T_E_obs = autocorr(nino3,60);
ACF_T_C_obs = autocorr(nino4,60);
ACF_T_CE_obs = autocorr(nino34,60);
%hw_store_a = mean(Hov_H_a(west_node-No*2,:));
%u_store_a = mean(Hov_u_a(west_node-No*2,:));
%ACF_h_W_model=autocorr(hw_store_a,60);
%ACF_u_model=autocorr(u_store_a,60);
%ACF_h_W_obs = autocorr(h_W_obs,60);
%ACF_u_obs = autocorr(u_obs,60);

save_spec = zeros(61, num_piece);
index1all = T3_store_a;
for j = 1:num_piece
    index1=index1all((j)*480+[1:480]);
    save_spec(:,j) = autocorr(index1,60);
end
mean_model= mean(save_spec');
std_model = std(save_spec');
low = mean_model-std_model;
high = mean_model+std_model;
x_axis = [0:1/12:5,5:-1/12:0];
y_axis = [low fliplr(high)];  
subplot('Position',[posx(6),posy(2),posw,posh])
%subplot(2,3,6)
hold on
hh0 = fill(x_axis, y_axis,'b','facealpha',facealpha);
hh0.EdgeColor = 'none';
hh1 = plot([0:60]/12, ACF_T_E_obs,'r','linewidth',2);
%hh2 = plot([0:60]/12, ACF_T_E_model,'b','linewidth',2);
hh2 = plot([0:60]/12, mean_model,'b','linewidth',2);
box on
set(gca,'fontsize',14)
title('Nino 3')
xlim([0,5])
xlabel('Year')
ylim([-0.6,1])
 

index1all = T34_store_a;
for j = 1:num_piece
    index1=index1all((j)*480+[1:480]);
    save_spec(:,j) = autocorr(index1,60);
end
mean_model= mean(save_spec');
std_model = std(save_spec');
low = mean_model-std_model;
high = mean_model+std_model;
x_axis = [0:1/12:5,5:-1/12:0];
y_axis = [low fliplr(high)]; 
%subplot(2,3,5)
subplot('Position',[posx(5),posy(2),posw,posh])
hold on
hh0 = fill(x_axis, y_axis,'b','facealpha',facealpha);
hh0.EdgeColor = 'none';
plot([0:60]/12, ACF_T_CE_obs,'r','linewidth',2)
%plot([0:60]/12, ACF_T_CE_model,'b','linewidth',2)
plot([0:60]/12, mean_model,'b','linewidth',2)
box on
set(gca,'fontsize',14)
title({['(d) ACFs'], ['Nino 3.4']})
xlim([0,5])
xlabel('Year')
ylim([-0.6,1])
 

index1all = T4_store_a;
for j = 1:num_piece
    index1=index1all((j)*480+[1:480]);
    save_spec(:,j) = autocorr(index1,60);
end
mean_model= mean(save_spec');
std_model = std(save_spec');
low = mean_model-std_model;
high = mean_model+std_model;
x_axis = [0:1/12:5,5:-1/12:0];
y_axis = [low fliplr(high)]; 
%subplot(2,3,4)
subplot('Position',[posx(4),posy(2),posw,posh])
hold on
hh0 = fill(x_axis, y_axis,'b','facealpha',facealpha);
hh0.EdgeColor = 'none';
plot([0:60]/12, ACF_T_C_obs,'r','linewidth',2)
%plot([0:60]/12, ACF_T_C_model,'b','linewidth',2)
plot([0:60]/12, mean_model,'b','linewidth',2)
box on
title('Nino 4')
xlim([0,5])
set(gca,'fontsize',14)
xlabel('Year')
ylim([-0.6,1])
ylabel('correlation')
 