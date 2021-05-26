clear

cte = set_cte();

for i = 1:length(cte.f)   
    %% Parameters
    f = cte.f(i);
    w = 2*pi*f;
    lambda = cte.c/f;
    k = w/cte.c; % 2*pi/lambda
    S = pi*(cte.D/2)^2;
    
    %% Side branch resonators
    % Helmholtz resonator:
    % resonator 1 (f = 40 Hz):
    f_HR = 40;
    [IL.helmholtz1(i),TL.helmholtz1(i), Z_HR_40(i)] = HR_IL_TL(f_HR, cte, i);
    TM(:,:,i) = [1,0; 1/Z_HR_40(i) 1/cte.rho_air];
    TL.HR(i) = 20*log10( abs(TM(1,1,i) + (S/cte.c)*TM(1,2,i) + (cte.c/S)*TM(2,1,i) + TM(2,2,i) )/2 );
    
    % resonator 1 (f = 1640 Hz):
    f_HR = 2000;
    [IL.helmholtz2(i),TL.helmholtz2(i), Z_HR_2000(i)] = HR_IL_TL(f_HR, cte, i);
    
    % Quarter-wavelength resonator:
    % 820Hz: 0.104mm
    % 1640Hz: 0.052mm
    H = 0.052; % lambda/4 [m]
    D_neck = 0.5*cte.D;
    
    S_1 = pi*(cte.D/2)^2;
    S_s = pi*(D_neck/2)^2;
    Z_s(i) = -1i*cte.rho_air*cte.c*cot(k*H);  % impedance lambda/4 resonator
    %Q = (cte.rho_air*cte.c/R_s)*sqrt();

    IL.lambda4(i) = 0;
    TL.lambda4(i) = 10*log10((tan(k*H)^2+4*(S_1/S_s)^2)/(4*(S_1/S_s)^2)); %-20*log(abs(2/(2+1i*(S_s/S_1)*tan(k*H)))); % same formula as before, but simplified
 
    %% TMM 
    p_i = 10;
    v_xi = 10;
    p_o = 1.5;
    v_xo = 1;
    
    syms T11 T12 T21 T22
    eq1 = T11*p_i+T12*v_xi == p_o;
    eq2 = T11*p_o+T12*v_xo == p_i;
    eq3 = T21*p_i+T22*v_xi == v_xo;
    eq4 = T21*p_o+T22*v_xo == v_xi;
    
    [A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [T11, T12, T21, T22]);
    T = linsolve(A,B);
    
    TM = [T(1), T(2); T(3), T(4)];

end

%% Plotting
figure(1),
subplot(3,1,1), hold on
plot(cte.f, abs(TL.helmholtz1)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")
plot(cte.f, abs(TL.helmholtz2)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")
subplot(3,1,2)
plot(cte.f, abs(TL.lambda4)), xlabel("f [Hz]"), ylabel("TR - lambda/4 [dB]")
subplot(3,1,3), hold on
plot(cte.f, abs(TL.helmholtz1)), xlabel("f [Hz]"),
plot(cte.f, abs(TL.helmholtz2)), xlabel("f [Hz]"),
plot(cte.f, abs(TL.lambda4)), xlabel("f [Hz]"), ylabel("TR [dB]")
legend("Helmholtz resonator 1","Helmholtz resonator 2", "lambda/4 ")

figure(3)
plot(cte.f, abs(TL.HR), cte.f, abs(TL.helmholtz1))


function [IL, TL, Z_HR] = HR_IL_TL(f_HR, cte, i)
    f = cte.f(i);
    w = 2*pi*f;
    k = w/cte.c; 
    
    D_neck = cte.D*0.2; 
    S_1 = pi*(cte.D/2)^2;
    S_s = pi*(D_neck/2)^2;
    
    l = (0.001:0.001:0.080);
    for j = 1:length(l)
        h = 0.080-l(j);
        V = (cte.c/(2*pi*f_HR))^2*(S_s/l(j));
        D_vol = sqrt((4*V)/(pi*h));
        if D_vol > 0.300
            D_vol = 0.300;
            h = (4/pi)*V/D_vol^2;
        end
        A_mat(j) = pi*D_neck*l(j) + pi*D_vol*h+2*(pi/4)*D_vol^2-(pi/4)*D_neck^2;
    end
    
    [A,I] = min(A_mat);
    l = l(I); % length where minimal material is used
    h = 0.080-l;
    
    V = (cte.c/(2*pi*f_HR))^2*(S_s/l);
    D_vol = sqrt((4*V)/(pi*h));
    S_vol = pi*(D_vol/2)^2;

    %Z_HR(i) = (1/S_s)*(1i*w*cte.rho_air*l*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator
    Z_HR = (1i*cte.rho_air*cte.c)*(S_vol*tan(k*l)*tan(k*h)-S_s)/(S_vol*tan(k*h)+S_s*tan(k*l)); % better impedance Helmholtz Resonator
    
    f_res = (cte.c/(2*pi))*(D_neck/D_vol)*sqrt(1/(l*h));
    
    IL = 0;
    TL = 20*log10(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR));
    
    % plot HZ-dimensions;
    if i == 1 && f_HR == 40
        disp(strcat('l: ',num2str(l)))
        disp(strcat('h: ',num2str(h)))
        disp(strcat('D_vol: ',num2str(D_vol)))
        
        figure(2),title("HR design"), hold on
        plot([-D_neck/2, -D_neck/2], [0 l],'b-', 'LineWidth', 2),
        plot([D_neck/2, D_neck/2], [0 l],'b-', 'LineWidth', 2),
        plot([-D_vol/2, -D_neck/2],[l l],'b-', 'LineWidth', 2),
        plot([D_neck/2 D_vol/2],[l l],'b-', 'LineWidth', 2),
        plot([-D_vol/2, -D_vol/2], [l l+h],'b-', 'LineWidth', 2),
        plot([D_vol/2, D_vol/2], [l l+h],'b-', 'LineWidth', 2),
        plot([-D_vol/2, D_vol/2],[l+h l+h],'b-', 'LineWidth', 2),
        axis([-0.3 0.3 -0.3 0.3])
        
    end
end
