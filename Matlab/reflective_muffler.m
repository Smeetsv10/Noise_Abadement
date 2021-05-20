function [IL,TL] = reflective_muffler(cte)
%% Reactive type Mufflers
    %% Import data
    % With expansion tube
    mic.A = read_table(readtable('mic1.csv','NumHeaderLines',1));
    mic.B = read_table(readtable('mic2.csv','NumHeaderLines',1));
    mic.C = read_table(readtable('mic3.csv','NumHeaderLines',1));
    mic.D = read_table(readtable('mic4.csv','NumHeaderLines',1));
    % Without expansion tube
    mic.A_wo = read_table(readtable('mic1n','NumHeaderLines',1));
    mic.B_wo = read_table(readtable('mic2n','NumHeaderLines',1));
    mic.C_wo = read_table(readtable('mic3n','NumHeaderLines',1));
    mic.D_wo = read_table(readtable('mic4n','NumHeaderLines',1));
    
for i = 1:length(cte.f)
    %% Parameters
    f = cte.f(i);
    w = 2*pi*f;
    lambda = cte.c/f;
    k = w/cte.c; % 2*pi/lambda
    
    %% Change in duct crosssection
    % Expansion chamber:
    D_new = cte.D*5; % upper limit is factor 5, 0.200m
    cte.L = 0.213;
    cte.m1 = 0.0785;
    cte.s1 = 0.010;
    cte.m2 = 0.0785;
    cte.s2 = 0.010;
    cte.d = 0.0935;
    
    A1 = (mic.A.p(i)*exp(1i*mic.A.phase(i))*exp(1i*k*cte.s1)-mic.B.p(i)*exp(1i*mic.B.phase(i)))*exp(-1i*k*(cte.m1+cte.s1))/(exp(1i*k*cte.s1)-exp(-1i*k*cte.s1)); % amp inlet
    A3 = (mic.C.p(i)*exp(1i*mic.C.phase(i))*exp(1i*k*cte.s2)-mic.D.p(i)*exp(1i*mic.D.phase(i)))*exp(-1i*k*(cte.m2-cte.d))/(exp(1i*k*cte.s2)-exp(-1i*k*cte.s2)); % amp outlet
    A3_wo = (mic.C_wo.p(i)*exp(1i*mic.C_wo.phase(i))*exp(1i*k*cte.s2)-mic.D_wo.p(i)*exp(1i*mic.D_wo.phase(i)))*exp(-1i*k*(cte.m2-cte.d))/(exp(1i*k*cte.s2)-exp(-1i*k*cte.s2));

    S_1 = pi*(cte.D/2)^2; % [m^2]
    S_2 = pi*(D_new/2)^2; % [m^2]
    N = S_1/S_2;

    R = 1; % reflection coefficient at the end of the duct
    IL.expansion(i) = 20*log10(abs( cos(k*cte.L) + 1i*0.5*(((1/N)+N)+R*((1/N)-N)*exp(-1i*2*k*cte.d))*sin(k*cte.L) )); % same formula as TL for R=0 ???
    IL.expansion_NX(i) = 20*log10(abs(A3_wo/A3));
    
    TL.expansion(i) = 10*log10(cos(k*cte.L)^2+0.25*(N+(1/N))^2*sin(k*cte.L)^2); % cos and + or 1 and - !
    TL.expansion_NX(i) = 10*log10(abs(A1/A3)^2); %10*log(abs((mic_A.p(i)/mic_C.p(i))^2))
        
    %% Side branch resonators
    % Helmholtz resonator:
    % resonator 1 (f = 1640 Hz):
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
    TL.helmholtz1(i) = 20*log10(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    
    % resonator 2 (f = 820 Hz):
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
    TL.helmholtz2(i) = 20*log10(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    
    % Quarter-wavelength resonator:
    % 820Hz: 0.104mm
    % 1640Hz: 0.052mm
    H = 0.052; % lambda/4 [m]
    D_neck = 0.5*cte.D;
    
    S_s = pi*(D_neck/2)^2;
    Z_s(i) = -1i*cte.rho_air*cte.c*cot(k*H);  % impedance lambda/4 resonator
    %Q = (cte.rho_air*cte.c/R_s)*sqrt();

    IL.lambda4(i) = 0;
    TL.lambda4(i) = 10*log10((tan(k*H)^2+4*(S_1/S_s)^2)/(4*(S_1/S_s)^2)); %-20*log(abs(2/(2+1i*(S_s/S_1)*tan(k*H)))); % same formula as before, but simplified

end

muffler_design(cte)

end

