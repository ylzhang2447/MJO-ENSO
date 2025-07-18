% Eigenvalues of the linearized skeleton model



% Skeleton model
Q_tilde = 0.9;
Gamma = 1.66*10;
epsilon = 0.0;
d_bar = 1;
H_bar = 22;
A_bar = 1;
c_o = 50;
L = 80/3/10;
dim_t = 34;
beta = 2.28*10^(-11);
Eig_Store = zeros(4,11);
Ev1 = zeros(4,11);
Ev2 = zeros(4,11);
Ev3 = zeros(4,11);
Ev4 = zeros(4,11);
eg_K_MJO = zeros(1,3);
eg_R_MJO = zeros(1,3);
eg_Q_MJO = zeros(1,3);
eg_A_MJO = zeros(1,3);
eg_K_Rossby = zeros(1,3);
eg_R_Rossby = zeros(1,3);
eg_Q_Rossby = zeros(1,3);
eg_A_Rossby = zeros(1,3);
for kk = 0:10
    if kk == 0
        kk = 0.01;
    end
    M_Omega = [1, 0, 0, 0;
        0, -1/3, 0, 0;        
        Q_tilde, -Q_tilde/3, 0, 0;
        0, 0, 0, 0];
    M_D = [-epsilon * d_bar, 0, 0, -H_bar/2;
        0, -epsilon * d_bar, 0, -H_bar/3;        
        0, 0, -epsilon * d_bar, (-1+Q_tilde/6) * H_bar;
        0, 0, Gamma * A_bar, 0];
    M = M_Omega * 1i * 10*kk/L*pi*2 - M_D;
    M = M;
    [VecEig, Lambda] = eig(M);
    if kk == 0.01
        kk = 0;
    end
    [sorting, ranking] = sort(diag(imag(Lambda)),'descend');
    Eig_Store(:,kk+1) = sorting;
    Ev1(:,kk+1) = sqrt((real(VecEig(:,ranking(1)))).^2 + (imag(VecEig(:,ranking(1)))).^2); %compute the norm of the eigenvectors
    Ev2(:,kk+1) = sqrt((real(VecEig(:,ranking(2)))).^2 + (imag(VecEig(:,ranking(2)))).^2);
    Ev3(:,kk+1) = sqrt((real(VecEig(:,ranking(3)))).^2 + (imag(VecEig(:,ranking(3)))).^2);
    Ev4(:,kk+1) = sqrt((real(VecEig(:,ranking(4)))).^2 + (imag(VecEig(:,ranking(4)))).^2);
    if kk >=1 && kk<=3
        eg_K_MJO(kk) = VecEig(1,ranking(2));
        eg_R_MJO(kk) = VecEig(2,ranking(2));
        eg_Q_MJO(kk) = VecEig(3,ranking(2));
        eg_A_MJO(kk) = VecEig(4,ranking(2));
        eg_K_Rossby(kk) = VecEig(1,ranking(3));
        eg_R_Rossby(kk) = VecEig(2,ranking(3));
        eg_Q_Rossby(kk) = VecEig(3,ranking(3));
        eg_A_Rossby(kk) = VecEig(4,ranking(3));
    end
end
lbd = diag(imag(Lambda));
lbd(ranking(2))
VecEig(:,ranking(2))

 
[xx,yy] = meshgrid(0:10,1:4);
fig1 = figure;
KK = [0.01,1:10];
KK1 = [1:10];
for i = 1:4
    subplot(4,3,3*i-1)
    plot(KK, Eig_Store(i,:)/2/pi/dim_t ,'bo-','linewidth',2);    
    set(gca,'fontsize',12)
    if i == 1
        title('Frequency (cpd)','fontsize',12)
    end
    if i == 4
        xlabel('Wavenumber','fontsize',12);
    end
    subplot(4,3,3*i-2);
    plot(KK1, Eig_Store(i,2:end)./(10*KK1/L*pi*2) * c_o,'bo-','linewidth',2);
    set(gca,'fontsize',12)
    if i == 1
        ylabel('Dry Kelvin','fontsize',12);
    elseif i == 2
        ylabel('MJO','fontsize',12);
    elseif i == 3
        ylabel('Moist Rossby','fontsize',12);
    elseif i == 4
        ylabel('Dry Rossby','fontsize',12);
    end
    if i == 1
        title('Phase speed (m/s)','fontsize',12)        
    end
    if i == 4
        xlabel('Wavenumber','fontsize',12);
    end
end
subplot(4,3,3);
amp = sqrt(real(Ev1).^2 + imag(Ev1).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
title('Eigenvectors','fontsize',12)
subplot(4,3,6);
amp = sqrt(real(Ev2).^2 + imag(Ev2).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
subplot(4,3,9);
amp = sqrt(real(Ev3).^2 + imag(Ev3).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
subplot(4,3,12);
amp = sqrt(real(Ev4).^2 + imag(Ev4).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
colormap jet

% saveas(fig1,"figures\old_MJO_eigen.png")
% saveas(fig1,"figures\old_MJO_eigen.eps",'epsc')