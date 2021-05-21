function [IL,TL, power] = calculations_muffler(cte)
    %% Import data
    [mic, power] = read_mics();
    
for i = 1:length(cte.f)
   %% Parameters
    f = cte.f(i);
    w = 2*pi*f;
    lambda = cte.c/f;
    k = w/cte.c; % 2*pi/lambda
    
%% 
% Reflective type Mufflers   
% _________________________

    %% Expansion chamber:
    D_new = cte.D*5; % upper limit is factor 5, 0.200m

    cte.A3_wo = (mic.wo.C.p(i)*exp(1i*mic.wo.C.phase(i))*exp(1i*k*cte.s2)-mic.wo.D.p(i)*exp(1i*mic.wo.D.phase(i)))*exp(-1i*k*(cte.m2-cte.d))/(exp(1i*k*cte.s2)-exp(-1i*k*cte.s2));
    %A4 = (mic.expCh.C.p(i)*exp(1i*mic.expCh.C.phase(i))-A3*exp(1i*k*-(cte.d-cte.m2))) / (exp(-1i*k*-(cte.d-cte.m2))); % zou 0 moeten zijn 
    
    S_1 = pi*(cte.D/2)^2; % [m^2]
    S_2 = pi*(D_new/2)^2; % [m^2]
    N = S_1/S_2;

    R = 1; % reflection coefficient at the end of the duct
    IL.expansion(i) = 20*log10(abs( cos(k*cte.L) + 1i*0.5*(((1/N)+N)+R*((1/N)-N)*exp(-1i*2*k*cte.d))*sin(k*cte.L) )); % same formula as TL for R=0 ???
    TL.expansion(i) = 10*log10(cos(k*cte.L)^2+0.25*(N+(1/N))^2*sin(k*cte.L)^2); % cos and + or 1 and - !
    [IL.expansion_NX(i), TL.expansion_NX(i)] = calc_IL_TL(mic.expCh, cte, i);
    
    
    %% Inlet/outlet extension
    [IL.inlet_outlet_NX(i),TL.inlet_outlet_NX(i)] = calc_IL_TL(mic.io, cte, i);
    
    %% Partitioning expansion chamber
    [IL.partitioning_NX(i),TL.partitioning_NX(i)] = calc_IL_TL(mic.prt, cte, i);

    %% Vibro acoustics
    power.expansion.abs(i) = sqrt(power.expansion.real(i)^2+power.expansion.imag(i)^2);
    power.normal.abs(i) = sqrt(power.normal.real(i)^2+power.normal.imag(i)^2);
    power.expansion.dB(i) = 10*log10(power.expansion.abs(i) / power.normal.abs(i));
 
    [IL.vibro_NX(i),TL.vibro_NX(i)] = calc_IL_TL(mic.vibro, cte, i);

    power.steel05.dB(i) = 10*log10( sqrt(power.steel05.real(i)^2+power.steel05.imag(i)^2) / power.normal.abs(i));
    power.steel15.dB(i) = 10*log10( sqrt(power.steel15.real(i)^2+power.steel15.imag(i)^2) / power.normal.abs(i));
    power.alum15.dB(i) = 10*log10( sqrt(power.alum15.real(i)^2+power.alum15.imag(i)^2) / power.normal.abs(i));
    
    t(i) = (0.5 + 0.0050*i)*10^-3; %0.5..1.5
    Vol_mat(i) = 2*(pi/4)*(D_new^2-cte.D^2)*t(i) + 2*pi*D_new*(cte.L-2*t(i))*t(i);
    tot_steel_price(i) = cte.price_steel*cte.rho_steel*Vol_mat(i);
    tot_titanium_price(i) = cte.price_titanium*cte.rho_titanium*Vol_mat(i);
%% 
% Absorption type Mufflers   
% _________________________

    %% Absorbing material
    %[IL.ab_NX(i),TL.ab_NX(i)] = calc_IL_TL(mic.ab, cte, i);

    %% Perforated ducts


%% 
% Rest
% _________________________

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

%     figure(12),
%     plot(t*1000, tot_steel_price, t*1000, tot_titanium_price), xlabel("thickness [mm]"), ylabel("Price [eur]")
%     legend('Steel', 'Titanium')

end

