new_path = '';

factor = 1;

rng(100)
% Parameters
La = 8/3;
LaDim = 40000;
LoDim = 17500;
dn = 2; % a factor that can refine the grids
Na = 64 * dn; % atmosphere grid points
No = round(LoDim/LaDim*Na);
Lo = LoDim/LaDim*La;
dx = La/Na; % distance between every two grid points
dim_t = 34*factor; % one dimension unit of time (0.34 days) 
x = (0:Na-1) * dx; % grid points in the x axis
xo = (0:No-1) * dx;
dt = 8/24/dim_t/10; % time step, every dt*dim_t hours  
gap = 30; % every gap points of saving the results
Gamma = 166*factor*0.1; % Gamma in skeleton model
Qbar = 0.9; % Q_bar in skeleton model
ee = 1;
cc1 = 0.18; % cc1 = 0.5;
rE = 1;
rW = 0.5;
Dim_T = 1.5;
Dim_H = 20.8;
Dim_U = 0.5;

zeta = 8.7*2.34;
gammaa = 6.529*0.9;
chia = 0.3086; % projection coefficient to atm
chio = 1.3801; % projection coefficient to ocean 

c = 0.05;
yy = -8:0.1:8;

% Atmopshere
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
 
Hbar = 22*factor;


aaa = 1;

wp = zeros(Na,1);
for i = 1:Na/2
    wp(i) = 0.88 +2.64*(1 ./ (1 + exp(-10*(x(i)+0.6))) - 1 ./ (1 + exp(-10*(x(i)-0.3))));
end
for i = Na/2+ 1:Na
    wp(i) = 0.88 +2.64*(1 ./ (1 + exp(-10*((x(i)-La)+0.6))) - 1 ./ (1 + exp(-10*((x(i)-La)-0.3))));
end

s_q = wp*aaa*factor;
s_theta = wp*aaa*factor;

S_q = s_q;
S_theta = s_theta;


Time_total = 360*10; % total integration time days 
N = round(Time_total/dim_t/dt); % total numerical integration time step

k = 2 * pi / (Na * dx) * (-Na/2:Na/2-1); % Fourier
%lk = 2 * pi * k / L;
% zero-th Fourier mode
if mod(Na,2) == 0
    zero_mode = Na/2+1;
else
    zero_mode = (Na+1)/2;
end

% multiplicative noise coefficient
lambda_A = 1.1*factor;
lambda_gsde = 1.1;
d_k = 1.1;
k_sde = 2;
multiplicative_factor = 2;
% initial values; cannot all be zero since otherwise the multiplicative
% noise will not be triggered
Ka = zeros(Na,N/gap);
Ra = zeros(Na,N/gap);
Ua = zeros(Na,N/gap);
Ka_bar = zeros(Na,N/gap);
Ra_bar = zeros(Na,N/gap);
tmp = sin((1/No:1/No:1)*2*pi);
T = zeros(No,N/gap); T(:,1) = 0.5*tmp;
Ko = zeros(No,N/gap); Ko(:,1) = 0.5*tmp;
Ro = zeros(No,N/gap); Ro(:,1) = -0.25*tmp;
K_phy = zeros(Na,N/gap);
R_phy = zeros(Na,N/gap);
Q = zeros(Na,N/gap); Q(Na,1) = 0.0;
Q_phy = zeros(Na,N/gap); Q_phy(:,1) = 0.1;
A = zeros(Na,N/gap); A(Na,1) = 0.001;
A_phy = zeros(Na,N/gap); A_phy(:,1) = 0.001;
A_bar_record = zeros(Na,N/gap);
Z = zeros(Na,N/gap);
Z_phy = zeros(Na,N/gap);
I_record = zeros(1,N/gap);
ap = zeros(1,N/gap);

noise_temp_record = zeros(Na,N/gap);

K_old = Ka(:,1);
R_old = Ra(:,1);
Q_old = Q(:,1);
A_old = A(:,1);
A_bar_old = A(:,1);
Z_old = Z(:,1);
K_phy_old = K_phy(:,1);
R_phy_old = R_phy(:,1);
Q_phy_old = Q_phy(:,1);
A_phy_old = A_phy(:,1);
Z_phy_old = Z_phy(:,1);
j = 1;
coeff = ones(1,N)*0.01;

west_node = 1:round(No*1/2);
east_node = round(No*1/2)+1:No;
T3_node = 31:53;
T4_node = 14:31;
T34_node = 25:43;


% Eq parameters
qc = 3.5; % latent heating multiplier coefficient
qe = 0.09296; % latent heating exponential coefficient
Tbar = 25/1.5; %1.5 is dimension, mean SST
tauq = 15; % latent heating adjustment rate
alpha_q = qc*qe*exp(qe*Tbar)/tauq;% latent heating factor

eta = (1.3 + 1.1 * tanh(7.5 * (xo - Lo/3)))*0.95; % Shifted east and weakened for shorter EP
eta2 = max(0, 4* exp(-((xo - Lo / (7/3)).^2 / 0.1)))*0.9; %


T_old = T(:,1);
A_bar = A(:,1);
Ko_old = Ko(:,1);
Ro_old = Ro(:,1);

% Prescribe the distribution of the decadal variable. Now we use a bimodal. 
xx = 0.01:0.01:1; 
p = ones(1,length(xx));
p = p/trapz(xx,p);

alpha = 0.5;
beta = 0.05;
gamma = 0.1;

m = trapz(xx,xx.*p); % mean value

lambda = dim_t/360/5; % damping time of the decadal variability, now 5 years
n = length(xx);
Phi = zeros(1,n);
sgm = zeros(1,n);
% build the stochastic process with multiplicative noise for the decadal
% variability
for i = 2:n
    Phi(i) = trapz(xx(1:i), (xx(1:i)-m) .* p(1:i));
    sgm(i) = real(sqrt(2/p(i) * (-lambda * Phi(i))));
end
tt1 = 1;
% decadal variability
I = zeros(1,N);I(1)=.5;
% integration in time
for yy = 1:200
    ii=1;

    sp = exp(-45*(xo-Lo/4).^2); % profile of the wind bursts spatial
    
    sp = transpose(sp);
    %ap = zeros(N,1);

    dt_noise = dt/10;ap_temp = 0; % a smaller time step for wind bursts to better simulate them
    iiap = 1;
    rd_save = randn(1,N*10);
    if yy > 1
        filename_last = [new_path 'main_test5new_Itest_',num2str(yy-1),'.mat'];
        %filename_last = ['E:\Projects\ENSO-MJO\codes_mjo_enso\GCM_model\record_u5_new_',num2str(yy-1),'.mat'];
        load(filename_last)
        j = 1;
        Ko_old = Ko_initial;
        Ro_old = Ro_initial;
        T_old = T_initial;
        K_phy_old = Ka_initial;
        R_phy_old = Ra_initial;
        Z_phy_old = Z_initial;
        A_phy_old = A_initial;
        A_bar_old = A_bar_initial;
        ap_temp = ap_initial;
        %u0 = u0_initial; ap_old = ap_initial;
        I = zeros(1,N);
        temp = round(I_initial*100) + 1;
        if temp<0
            temp = 1;
        elseif temp > n
            temp = n;
        end
        sgm_x = sgm(temp);
        I(1) = I_initial + (-lambda * (I_initial - m) * dt) + sgm_x * randn *sqrt(dt);

        %initial
        Ua = zeros(Na,N/gap);
        T = zeros(No,N/gap);
        Ko = zeros(No,N/gap);
        Ro = zeros(No,N/gap);
        K_phy = zeros(Na,N/gap);
        R_phy = zeros(Na,N/gap);
        K = zeros(Na,N/gap);
        R = zeros(Na,N/gap);        
        Q = zeros(Na,N/gap);
        Q_phy = zeros(Na,N/gap);
        A = zeros(Na,N/gap);
        A_phy = zeros(Na,N/gap);
        A_bar_record = zeros(Na,N/gap);
        ap = zeros(1,N/gap);
    end


    for i = 2:N
        temp = round(I(i-1)*100) + 1;
        if temp<0
            temp = 1;
        elseif temp > n
            temp = n;
        end
        sgm_x = sgm(temp);
        I(i) = I(i-1) + (-lambda * (I(i-1) - m) * dt) + sgm_x * randn *sqrt(dt);
    end
    I = 0.25 + I*0.5;

    for i = 1:N

        Tw = mean(T_old(west_node));
        Te = mean(T_old(east_node));
        T4 = mean(T_old(T4_node));
        T3 = mean(T_old(T3_node));
        phase = 0;

        temp_alphaq = (1.8 + 1*(- eta2/3 + (0.2 + abs(mean(T_old(T4_node))+0.4).* eta2 ).^2/5));
        
        season_alphaq = 1. + (0.4*sin(2*pi*(dt*tt1*dim_t/360+phase)) +  0.5*sin(2*pi*(dt*tt1*dim_t/360+phase-2/12)) * eta2/3.6 + 0.2*sin(2*pi*(dt*tt1*dim_t/360+1/12))*eta/2.4/0.95);
        
   
        tt1 = tt1+1;

        Eq = zeros(Na,1); 
        Eq2 = zeros(Na,1); 

        Eq(1:No) = alpha_q*(T_old);

        Eq2(1:No) = 20*alpha_q * (2.3 + 0.7*eta/0.95)' .* T_old;

        Eqzma = Eq-mean(Eq);
        Eqzma2 = Eq2-mean(Eq2);

        A_bar = 10 / Hbar / (1 - Qbar) * (Eqzma*chia + S_q - Qbar * S_theta); % mean state of A; in coupled ENSO-MJO, feedback from Eq should be added
        
        temp_H = -(3/2)*(Eqzma)*chia/(1-Qbar);

        da = 10e-8;

        AA=zeros(Na,Na);
        for ix=1:Na-1 
            AA(ix,ix+1)= 1/dx; 
        end
        AA(Na,1)=1/dx; % A(i,i+1)
        for ix=1:Na
            AA(ix,ix)= da -1/dx; 
        end
        % Compute vector B
        BB1=zeros(Na,1);
        for ix=1:Na
            BB1(ix,1)= temp_H(ix,1); 
        end 
        BB1(1:Na,1)=BB1(1:Na,1)-mean(BB1(1:Na,1));% remove zonal average
        % Compute X (matrix inversion) 
        XX=AA\BB1;


        T_temp = zeros(Na,1); T_temp(1:No) = T_old;

        temp_wind = 0;

        force = -1/2 * (Hbar * A_phy_old - temp_wind);
        c_k = 10;
        K_old = fftshift(fft(K_phy_old));
        force = fftshift(fft(force)); 
        
        K_new = K_old .* exp( - d_k * dt - 1i * k' * c_k * dt) + ( force ./ ( d_k + 1i * k' * c_k) ) .* ( 1 - exp( - d_k * dt - 1i * k' * c_k * dt ) );
        if d_k == 0
            K_new(zero_mode) = K_old(zero_mode) + force(zero_mode) * dt;
        else
            K_new(zero_mode) = K_old(zero_mode) * exp( - d_k * dt ) + force(zero_mode) / d_k * ( 1 - exp( - d_k * dt ) );
        end
        K_phy_new = ifft(ifftshift(K_new));
        K_phy_new = real(K_phy_new); 
        
        
        % update R in Fourier space
        force = -1/3 * (Hbar * A_phy_old + temp_wind);
        c_k = -10/3;
        R_old = fftshift(fft(R_phy_old));
        force = fftshift(fft(force)); 
        
        R_new = R_old .* exp( - d_k * dt - 1i * k' * c_k * dt) + ( force ./ ( d_k + 1i * k' * c_k) ) .* ( 1 - exp( - d_k * dt - 1i * k' * c_k * dt ) );
        if d_k == 0
            R_new(zero_mode) = R_old(zero_mode) + force(zero_mode) * dt;
        else
            R_new(zero_mode) = R_old(zero_mode) * exp( - d_k * dt ) + force(zero_mode) / d_k * ( 1 - exp( - d_k * dt ) );
        end
        R_phy_new = ifft(ifftshift(R_new));
        R_phy_new = real(R_phy_new); 

        season_temp = 1;
       
        sp_temp = zeros(size(Eq));
        sp_temp(1:No) = sp;
        %temp_Eq = Eqzma + 0.1*S_q*(1-I(i));
        temp_Eq = (Eqzma2) + (S_q)*(1-1*I(i));

        noise_Q = max(temp_Eq, 0.1);

        noise_Q =  4*(0.77+tanh(noise_Q*2.3*0.1-1))*season_temp;
        noise_temp = noise_Q;



        Z_phy_new = Z_phy_old + ( - d_k * Z_phy_old - Hbar * (1 - Qbar) * A_phy_old ) * dt + sqrt(dt)*noise_temp.*randn(Na,1);
        Q_phy_new = Z_phy_new + Qbar * ( K_phy_new + R_phy_new);
        Q_new = fftshift(fft(Q_phy_new));

        % update A in physical space
        A_total = A_phy_old + A_bar; % total convective activity A
        A_total(A_total <= 0) = 10e-5; % positivity
        term1_A = Gamma * A_total .* Q_phy_old; % term 1: dynamical part
        temp_q = Q_phy_old; temp_q(temp_q<0) = 0;
        temp_a = A_bar; temp_a(temp_a<0) = 0;
        term2_A = - lambda_A*A_phy_old; % term 2: damping
        term3_A =  sqrt(10*Gamma*abs(Q_phy_old).*A_total).*randn(Na,1);
        A_temp = A_phy_old + (term1_A + term2_A) * dt + term3_A * sqrt(dt);
        A_temp(A_temp <= -A_bar) = - A_bar(A_temp <= -A_bar) + 10e-5; % positivity
        A_phy_new = A_temp;
        A_new = fftshift(fft(A_phy_new));
        


        taux =  gammaa*(XX + (K_phy_old - R_phy_old));

        K_old = K_new; K_phy_old = K_phy_new; 
        R_old = R_new; R_phy_old = R_phy_new; 
        Q_old = Q_new; Q_phy_old = Q_phy_new; 
        A_old = A_new; A_phy_old = A_phy_new; 
        Z_phy_old = Z_phy_new; 

        
        % approximation of a cubic damping together with the seasonal effect
        b = cc1 * zeta * alpha_q *temp_alphaq.*season_alphaq;
        

        tempK = zeros(No,1);
        tempR = zeros(No,1);
        
        for ix=2:No
            tempK(ix,1) = - ee * cc1/dx * (Ko_old(ix,1)-Ko_old(ix-1,1)) + ee * chio*cc1/2 * taux(ix,1); 
        end
        tempK(1,1) = - ee*cc1/dx*(Ko_old(1,1)-rW*Ro_old(1,1)) + ee*chio*cc1/2*taux(1,1);
        for ix=1:No-1
            tempR(ix,1) = ee*cc1/(3*dx)*(Ro_old(ix+1,1)-Ro_old(ix,1)) - ee*chio*cc1/3*taux(ix,1); 
        end
        tempR(No,1) = ee*cc1/(3*dx)*(rE*Ko_old(No,1)-Ro_old(No,1)) - ee*chio*cc1/3*taux(No,1);
        
        Ko_new=Ko_old+dt*tempK;
        Ro_new=Ro_old+dt*tempR;
        
        c2 = 0.1;
        tempT= -ee*b'.*T_old + ee*cc1*eta'.*(Ko_old +Ro_old) + (0.1 + 1.2*I(i))*ee*cc1*eta2'.*(Ko_old-Ro_old) + 1*ee*cc1*(1.5*eta2'-1.*eta'/0.95)* chio*c2;
        T_new = T_old + dt*tempT;

    

        T_old = T_new;
        Ko_old = Ko_new;
        Ro_old = Ro_new;
        A_bar_old = A_bar;
        
        %A_new = A_bar;    
        if mod(i,gap) == 0
            Ua(:,j) = XX;
            %Ka_bar(:,j) = XX2;
            %Ra_bar(:,j) = XX3;
            A(:,j) = A_new;
            Ko(:,j) = Ko_new;
            Ro(:,j) = Ro_new;
            T(:,j) = T_new;
            I_record(j) = I(i);
            ap(j) = ap_temp;
            %
            K(:,j) = K_new;
            R(:,j) = R_new;
            Q(:,j) = Q_new;
            K_phy(:,j) = K_phy_new;
            R_phy(:,j) = R_phy_new;
            Q_phy(:,j) = Q_phy_new;
            A_phy(:,j) = A_phy_new;
            A_bar_record(:,j) = A_bar;

            noise_temp_record(:,j) = noise_temp;
            %Eq_record(:,j) = temp_Eq;
            if j == N/gap
               Ko_initial = Ko_new;
               Ro_initial = Ro_new;
               T_initial = T_new;
               I_initial = I(i);
               Ka_initial = K_phy_new;
               Ra_initial = R_phy_new;
               Z_initial = Z_phy_new;
               A_initial = A_phy_new;
               A_bar_initial = A_bar;
               ap_initial = ap_temp;

               filenames = [new_path 'main_test5new_Itest_',num2str(yy),'.mat'];
               %filenames = ['E:\Projects\ENSO-MJO\codes_mjo_enso\GCM_model\record_u5_new_',num2str(yy),'.mat'];
               disp(filenames)
               %save(filenames,'Eq_record','noise_temp_record','Ua','Ka_bar','Ra_bar', 'A','Ko','Ro','T','I_record','A_bar_record','Q_phy','A_phy','R_phy','K_phy','A_bar_initial','A_initial','Ko_initial','Ro_initial','T_initial','Ka_initial','Z_initial','Ra_initial','I_initial')
               save(filenames,'ap','ap_initial','noise_temp_record','Ua','A','Ko','Ro','T','I_record','A_bar_record','Q_phy','A_phy','R_phy','K_phy','A_bar_initial','A_initial','Ko_initial','Ro_initial','T_initial','Ka_initial','Z_initial','Ra_initial','I_initial')
            end
            j = j + 1;
        end
    end
end
