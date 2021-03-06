clear
close all

cte = set_cte();
mic = read_mics();

% 40 Hz
vel.helm1.B = read_vel(readtable('vel2h.csv','NumHeaderLines',1));
vel.helm1.C = read_vel(readtable('vel3h.csv','NumHeaderLines',1));
vel.helmv1.B = read_vel(readtable('vel2hv.csv','NumHeaderLines',1));
vel.helmv1.C = read_vel(readtable('vel3hv.csv','NumHeaderLines',1));

vel.helm_new.B = read_vel(readtable('velbh.csv','NumHeaderLines',1));
vel.helm_new.C = read_vel(readtable('velch.csv','NumHeaderLines',1));
vel.helm_newv.B = read_vel(readtable('velbhv.csv','NumHeaderLines',1));
vel.helm_newv.C = read_vel(readtable('velchv.csv','NumHeaderLines',1));
% array close
vel.helm2.B = read_vel(readtable('vel2h1.csv','NumHeaderLines',1));
vel.helm2.C = read_vel(readtable('vel3h1.csv','NumHeaderLines',1));
vel.helmv2.B = read_vel(readtable('vel2hv1.csv','NumHeaderLines',1));
vel.helmv2.C = read_vel(readtable('vel3hv1.csv','NumHeaderLines',1));
% array further
vel.helm3.B = read_vel(readtable('vel2h2.csv','NumHeaderLines',1));
vel.helm3.C = read_vel(readtable('vel3h2.csv','NumHeaderLines',1));
vel.helmv3.B = read_vel(readtable('vel2hv2.csv','NumHeaderLines',1));
vel.helmv3.C = read_vel(readtable('vel3hv2.csv','NumHeaderLines',1));
% 1 parallel 40 and 50 Hz
vel.helm4.B = read_vel(readtable('vel2pa.csv','NumHeaderLines',1));
vel.helm4.C = read_vel(readtable('vel3pa.csv','NumHeaderLines',1));
vel.helmv4.B = read_vel(readtable('vel2pav.csv','NumHeaderLines',1));
vel.helmv4.C = read_vel(readtable('vel3pav.csv','NumHeaderLines',1));
% array parallel 40 and 50 Hz
vel.helm5.B = read_vel(readtable('vel2pa1.csv','NumHeaderLines',1));
vel.helm5.C = read_vel(readtable('vel3pa1.csv','NumHeaderLines',1));
vel.helmv5.B = read_vel(readtable('vel2pav1.csv','NumHeaderLines',1));
vel.helmv5.C = read_vel(readtable('vel3pav1.csv','NumHeaderLines',1));
%%

for i = 1:length(cte.fh)   
    %% Parameters
    f = cte.fh(i);
    w = 2*pi*f;
    lambda = cte.c/f;
    k = w/cte.c; % 2*pi/lambda
    
    S_1 = pi*(cte.D/2)^2;
    D_neck = cte.D*0.15; 
    S_s = pi*(D_neck/2)^2;
    f_HR = 40;
    
    %% 1 HR
%     l_bot = 0.038;
%     l = (l_bot:0.001:0.060); 
%     for j = 1:length(l)
%         h(j) = 0.062-l(j);
%         l_eq(j) = l(j) +1.7*D_neck/2;
%         V(j) = (cte.c/(2*pi*f_HR))^2*(S_s/l_eq(j)); % volume for f_HR
%         %V = width*(0.200*h + 2*0.078*(0.038+l));
%         width(j) = V(j)/(0.200*h(j)+2*0.078*(0.38+l(j)));
%     end

    l = 0.050;
    h = 0.060-l;
    l_eq = l +1.7*(D_neck/2);
    V = (cte.c/(2*pi*f_HR))^2*(S_s/l_eq); % volume for f_HR
    width = V/(0.200*h+2*0.078*(0.038+l_eq));
   
    R = 0.00;
    %Z_HR(i) = R + (1/S_s)*(1i*w*cte.rho_air*l_eq*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator
    Z_HR(i) = (1i*cte.rho_air)*(w*l_eq - (cte.c^2*S_s)/(w*V)); % impedance Helmholtz Resonator
    %Z_HR = (1i*cte.rho_air*cte.c)*(S_vol*tan(k*l)*tan(k*h)-S_s)/(S_vol*tan(k*h)+S_s*tan(k*l)); % better impedance Helmholtz Resonator
    
    f_res = (cte.c/(2*pi))*sqrt(S_s/(l_eq*V));

    TL.helmholtz1(i) = 20*log10(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    TM.helmholtz1(i,:,:) = [1, 0; 1/Z_HR(i) 1];
    
    %% Series HR
    M = cte.rho_air*S_s*l_eq;
    K = cte.rho_air*cte.c^2*S_s^2/V;
    a11 = -w^2*M+1i*w*R+K;
    a12 = -K;
    a21 = -K;
    a22 = -w^2*M+1i*w*R+ K;
    
    Z_HR2(i) = 1/(1i*w*S_s^2)*(a11*a22-a12*a21)/a22;
    TL.helmholtz2(i) = 20*log10(0.5*abs(2+(cte.rho_air*cte.c*(S_s/S_1))*(1/Z_HR2(i))));
    
    %% TMM
    % 40 Hz
    [TL.HR_40(i), TM.HR_40(i,:,:)] = TL_HR(mic.helm1, vel.helm1, i);
    %[TL.HR_40(i), TM.HR_40(i,:,:)] = TL_HR(mic.helm_new, vel.helm_new, i);
    
    cte.m1 = 0.065;
    cte.s1 = 0.010;
    cte.m2 = 0.320;
    cte.s2 = 0.010;
    cte.A3_wo = (mic.wo.C.p(i)*exp(1i*mic.wo.C.phase(i))*exp(1i*k*cte.s2)-mic.wo.D.p(i)*exp(1i*mic.wo.D.phase(i)))*exp(-1i*k*(cte.m2-cte.d))/(exp(1i*k*cte.s2)-exp(-1i*k*cte.s2));
    TL.HR_method(i) = TL_helmholtz(mic.helm1, cte, i);
    
    TM.HR(:,:) = TM.HR_40(i,:,:);
    TM.HR2(i,:,:) = TM.HR(:,:)^2;
    TL.HR2(i)= 20*log10(abs(  TM.HR2(i,1,1) + (S_s/cte.c)*TM.HR2(i,1,2)  + (cte.c/S_s)*TM.HR2(i,2,1) + TM.HR2(i,2,2) )/2);
    
    % array close
    [TL.HR_array_close(i), TM.HR_array_close(i,:,:)] = TL_HR(mic.helm2, vel.helm2, i);
    
    % array further
	[TL.HR_array_further(i), TM.HR_array_further(i,:,:)] = TL_HR(mic.helm3, vel.helm3, i); 
    
    % 1 parallel
	[TL.HR_parallel(i), TM.HR_parallel(i,:,:)] = TL_HR(mic.helm4, vel.helm4, i);
    
    % array parallel
	[TL.HR_parallel_array(i), TM.HR_parallel_array(i,:,:)] = TL_HR(mic.helm5, vel.helm5, i);
    
    % with phase
    pi1(i) = mic.helm_new.B.p(i)*10^6*exp(1i* mic.helm1.B.phase(i));
    po1(i) = mic.helm_new.C.p(i)*10^6*exp(1i* mic.helm1.C.phase(i));
    pi2(i) = mic.helm_newv.B.p(i)*10^6*exp(1i* mic.helmv1.B.phase(i));
    po2(i) = mic.helm_newv.C.p(i)*10^6*exp(1i* mic.helmv1.C.phase(i));
    vi1(i) = vel.helm_new.B.v(i)*exp(1i* vel.helm1.B.phase(i));
    vo1(i) = vel.helm_new.C.v(i)*exp(1i* vel.helm1.C.phase(i));
    vi2(i) = vel.helm_newv.B.v(i)*exp(1i* vel.helmv1.B.phase(i));
    vo2(i) = vel.helm_newv.C.v(i)*exp(1i* vel.helmv1.C.phase(i));
    P = [po1(i), po2(i); vo1(i), vo2(i)];
    R = [pi1(i), pi2(i); vi1(i), vi2(i)];
    T(i,:,:) = P/R;
    TL_withphase(i)= 20*log10(abs(  T(i,1,1) + (S_s/cte.c)*T(i,1,2)  + (cte.c/S_s)*T(i,2,1) + T(i,2,2) )/2);
    
    
    p1(i) = mic.helm_new.B.p(i)*10^6*exp(1i* mic.helm1.B.phase(i));
    p2(i) = mic.helm_new.C.p(i)*10^6*exp(1i* mic.helm1.C.phase(i));
    v1(i) = vel.helm_new.B.v(i)*exp(1i* vel.helm1.B.phase(i));
    v2(i) = vel.helm_new.C.v(i)*exp(1i* vel.helm1.C.phase(i));
    P = [p2(i), p1(i); v2(i), -v1(i)];
    R = [p1(i), p2(i); v1(i), -v2(i)];
    T(i,:,:) = P/R;
    TL_withphase(i)= 20*log10(abs(  T(i,1,1) + (S_s/cte.c)*T(i,1,2)  + (cte.c/S_s)*T(i,2,1) + T(i,2,2) )/2);
    
%     TM.HR(:,:) = T(i,:,:);
%     TM.HR2(i,:,:) = TM.HR(:,:)^2;
%     TL.HR2(i)= 20*log10(abs(  TM.HR2(i,1,1) + (S_s/cte.c)*TM.HR2(i,1,2)  + (cte.c/S_s)*TM.HR2(i,2,1) + TM.HR2(i,2,2) )/2);
    
    %[p1(i), p2(i)] = TL_helmholtz(mic.helm1, cte, i);
end

%% Plotting
disp(strcat('D_neck: ',num2str(D_neck)))
disp(strcat('l: ',num2str(l)))
disp(strcat('h: ',num2str(h)))
disp(strcat('width: ',num2str(width)))

figure(2), hold on
plot(cte.fh, abs(TL.HR_method))
plot(cte.fh, abs(TL_withphase)), xlabel("f [Hz]"), ylabel("TL - Helmholtz resonator [dB]")
legend('Analytical', 'Numerical')

% figure(3), hold on
% plot(cte.fh, TM.HR_40(:,1,2))
% plot(cte.fh, TM.HR_array_close(:,1,2))
% plot(cte.fh, TM.HR_array_further(:,1,2))
% xlabel("f [Hz]"), ylabel("T11")
% legend('40 Hz', 'Array close', 'Array further')

figure(4), hold on
plot(cte.fh, abs(TL.HR2))
plot(cte.fh, abs(TL.HR_array_close))
plot(cte.fh, abs(TL.HR_array_further))
xlabel("f [Hz]"), ylabel("TL - Array HR [dB]")
legend('40 Hz', 'Array close', 'Array further')

figure(5), hold on
plot(cte.fh, abs(TL.HR_parallel))
plot(cte.fh, abs(TL.HR_parallel_array))
xlabel("f [Hz]"), ylabel("TL - Parallel HR [dB]")
legend('1 parallel', 'Array parallel')

figure(6), hold on
plot(cte.fh, abs(TL_withphase))
plot(cte.fh, abs(TL.HR_40))
legend('mag and phase', 'mag only')


% figure(7)
% subplot(4,1,1), hold on
% plot(cte.fh, TM.HR2(:,1,1)),
% plot(cte.fh, TM.HR_array_close(:,1,1)),
% plot(cte.fh, TM.HR_array_further(:,1,1)),
% legend('40 Hz^2', 'Array close', 'Array further')
% subplot(4,1,2), hold on
% plot(cte.fh, TM.HR2(:,1,2)),
% plot(cte.fh, TM.HR_array_close(:,1,2)),
% plot(cte.fh, TM.HR_array_further(:,1,2)),
% legend('40 Hz^2', 'Array close', 'Array further')
% subplot(4,1,3), hold on
% plot(cte.fh, TM.HR2(:,2,1)),
% plot(cte.fh, TM.HR_array_close(:,2,1)),
% plot(cte.fh, TM.HR_array_further(:,2,1)),
% legend('40 Hz^2', 'Array close', 'Array further')
% subplot(4,1,4), hold on
% plot(cte.fh, TM.HR2(:,2,2))
% plot(cte.fh, TM.HR_array_close(:,2,2))
% plot(cte.fh, TM.HR_array_further(:,2,2))
% legend('40 Hz^2', 'Array close', 'Array further')

% figure(8), hold on
% plot(cte.fh, vel.helm1.B.v)
% plot(cte.fh, vel.helm_new.B.v)
% plot(cte.fh, vel.helm_newv.B.v)
% xlabel('f'), ylabel('velocity')
% legend('old', 'new', 'new v')

tilefigs

%% Functions
function table = read_vel(T)
    cte = set_cte();
    S_s = (pi/4)*cte.D^2;
    
    table.freq = T.x_Frequency_Hz_;
    table.v = cte.rho_air*S_s*T.Velocity_mm_sec__Magnitude*10^-3;
    table.phase = T.Angle_degrees__Phase*(pi/180); 
end

function [TL, TM] = TL_HR(mic, vel, i)
    cte = set_cte();
    f = cte.fh(i);
    w = 2*pi*f;
    D_neck = cte.D*0.15; 
    S_s = pi*(D_neck/2)^2;

    p1(i) = mic.B.p(i)*10^6;
    p2(i) = mic.C.p(i)*10^6;
    v1(i) = vel.B.v(i);
    v2(i) = vel.C.v(i);
    
    P = [p2(i), p1(i); v2(i), -v1(i)];
    R = [p1(i), p2(i); v1(i), -v2(i)];
    TM = P/R; % P = T*R
    TM_11_mag(i) = (TM(1,1));
    TM_12_mag(i) = (TM(1,2));
    TM_21_mag(i) = (TM(2,1));
    TM_22_mag(i) = (TM(2,2));
    
    TL = 20*log10(abs( TM_11_mag(i) + (S_s/cte.c)*TM_12_mag(i) + (cte.c/S_s)*TM_21_mag(i) + TM_22_mag(i))/2);
end
