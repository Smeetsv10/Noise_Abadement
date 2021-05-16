function [IL,TL] = reflective_muffler(cte)
%% Reactive type Mufflers
    %% Import data
    % With expansion tube
    mic_A = read_table(readtable('mic_5.csv','NumHeaderLines',1));
    mic_B = read_table(readtable('mic_6.csv','NumHeaderLines',1));
    mic_C = read_table(readtable('mic_7.csv','NumHeaderLines',1));
    mic_D = read_table(readtable('mic_8.csv','NumHeaderLines',1));
    % Without expansion tube
    mic_C_without = read_table(readtable('mic_C_without.csv','NumHeaderLines',1));
    
for i = 1:length(cte.f)
    %% Parameters
    f = cte.f(i);
    w = 2*pi*f;
    lambda = cte.c/f;
    k = w/cte.c; % 2*pi/lambda
    
    %% Change in duct crosssection
    % Expansion chamber:
    D_new = cte.D*5; % upper limit is factor 5, 0.200m
    L = 0.213;
    m1 = 0.122;
    s1 = 0.010;
    m2 = 0.035;
    s2 = 0.010;
    d = 0.050;
    
    A1 = (mic_A.p(i)*exp(1i*k*s1)-mic_B.p(i))*exp(-1i*k*(m1+s1))/(exp(1i*k*s1)-exp(-1i*k*s1)); % amp inlet
    A3 = (mic_C.p(i)*exp(1i*k*s2)-mic_D.p(i))*exp(-1i*k*(m2-d))/(exp(1i*k*s2)-exp(-1i*k*s2)); % amp outlet
    
    S_1 = pi*(cte.D/2)^2; % [m^2]
    S_2 = pi*(D_new/2)^2; % [m^2]
    N = S_1/S_2;
    
    %IL.expansion(i) = 20*log( abs(cos(k*L) + 1i*0.5*((1/N)+N+R*((1/N)-N)*exp(-1i*2*k*d))*sin(k*L)) );
    IL.expansion_NX(i) = 20*log(abs(mic_C_without.p(i)/mic_C.p(i)));
    TL.expansion(i) = 10*log(1+0.25*(N-(1/N))^2*sin(k*L)^2); % maybe replace 1 by cos(k*L)^2
    TL.expansion_NX(i) = 10*log(abs(A1/A3)^2); %10*log(abs((mic_A.p(i)/mic_C.p(i))^2))
        
    %% Side branch resonators
    % Helmholtz resonator:
    % resonator 1:
    l = 0.010; % m
    D_neck = cte.D*0.3;
    h = 0.070; % max 70
    D_vol = 0.030;
    V = (pi/4)*D_vol^2*h; % m^3

    S_1 = pi*(cte.D/2)^2;
    S_s = pi*(D_neck/2)^2;

    Z_HR(i) = (1/S_s)*(1i*w*cte.rho_air*l*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator

    f_res(i) = cte.c/(2*pi)*sqrt(S_s/(l*V));

    IL.helmholtz1(i) = 0;
    TL.helmholtz1(i) = 20*log(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    
    % resonator 2:
    l = 0.010; % m
    D_neck = cte.D*0.3;
    h = 0.070; % max 70
    D_vol = 0.015;
    V = (pi/4)*D_vol^2*h; % m^3

    S_1 = pi*(cte.D/2)^2;
    S_s = pi*(D_neck/2)^2;

    Z_HR(i) = (1/S_s)*(1i*w*cte.rho_air*l*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator

    f_res(i) = cte.c/(2*pi)*sqrt(S_s/(l*V));

    IL.helmholtz2(i) = 0;
    TL.helmholtz2(i) = 20*log(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    
    % Quarter-wavelength resonator:
    % 820Hz: 0.104mm
    % 1640Hz: 0.052mm
    H = 0.104; % lambda/4 [m]
    Z_s(i) = -1i*cte.rho_air*cte.c*cot(k*H);  % impedance lambda/4 resonator
    %Q = (cte.rho_air*cte.c/R_s)*sqrt();

    IL.lambda4(i) = 0;
    TL.lambda4(i) = -20*log(abs(2/(2+1i*(S_s/S_1)*tan(k*H)))); % same formula as before, but simplified

end

%% Post Processing
%     % Helmholtz resonator
%     figure(1),
%     plot(cte.f/transpose(f_res), abs(TL.helmholtz)), xlabel("f/f_{res}"), ylabel("TR - Helmholtz resonator")
%     
%     % Quarter wavelength resonator
%     figure(2),
%     plot(cte.f*H/cte.c, abs(TL.lambda4)), xlabel("H/lambda"), ylabel("TR - lambda/4")
% 
%     % Expansion chamber
%     figure(3),
%     plot(2*pi*cte.f/cte.c*L, abs(TL.expansion)), xlabel("kL"), ylabel("TR - expansion chamber")

end

