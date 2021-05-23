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
    S_vol = pi*(D_vol/2)^2;

    Z_HR(i) = (1/S_s)*(1i*w*cte.rho_air*l*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator
    Z_HR(i) = (1i*cte.rho_air*cte.c)*(S_vol*tan(k*l)*tan(k*h)-S_s)/(S_vol*tan(k*h)+S_s*tan(k*l)); % better impedance Helmholtz Resonator
    
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
 

if figures(1)
    
    figure(1),
    subplot(4,1,1), hold on
    plot(cte.f, abs(TL.helmholtz1)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")
    plot(cte.f, abs(TL.helmholtz2)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")
    subplot(4,1,2)
    plot(cte.f, abs(TL.lambda4)), xlabel("f [Hz]"), ylabel("TR - lambda/4 [dB]")
    subplot(4,1,3)
    plot(cte.f, abs(TL.expansion)), xlabel("f [Hz]"), ylabel("TR - expansion chamber [dB]")
    subplot(4,1,4), hold on
    plot(cte.f, abs(TL.helmholtz1)), xlabel("f [Hz]"),
    plot(cte.f, abs(TL.helmholtz2)), xlabel("f [Hz]"),
    plot(cte.f, abs(TL.expansion)), xlabel("f [Hz]"),
    plot(cte.f, abs(TL.lambda4)), xlabel("f [Hz]"), ylabel("TR [dB]")
    legend("Helmholtz resonator 1","Helmholtz resonator 2", "expansion chamber", "lambda/4 ")
        
    figure(2),
    TL.total = TL.expansion + TL.lambda4;
    plot(cte.f, abs(TL.total)), xlabel("f [Hz]"), ylabel("total TR [dB]")
end