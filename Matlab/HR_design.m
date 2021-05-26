clear

cte = set_cte();


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
    width = V/(0.200*h+2*0.078*(0.038+l));
   
    Z_HR(i) = 0.2 + (1/S_s)*(1i*w*cte.rho_air*l_eq*S_s + (cte.rho_air*cte.c^2*S_s^2)/(1i*w*V)); % impedance Helmholtz Resonator
    %Z_HR(i) = (1i*cte.rho_air)*(w*l - (cte.c^2*S_s)/(w*V)); % impedance Helmholtz Resonator
    %Z_HR = (1i*cte.rho_air*cte.c)*(S_vol*tan(k*l)*tan(k*h)-S_s)/(S_vol*tan(k*h)+S_s*tan(k*l)); % better impedance Helmholtz Resonator
    
    f_res = (cte.c/(2*pi))*sqrt(S_s/(l_eq*V));

    TL.helmholtz1(i) = 20*log10(abs(1+0.5*(S_s/S_1)*cte.rho_air*cte.c/Z_HR(i)));
    
    
    
   %% Plotting
   if i==1

   end
   
end

disp(strcat('D_neck: ',num2str(D_neck)))
disp(strcat('l: ',num2str(l)))
disp(strcat('h: ',num2str(h)))
disp(strcat('width: ',num2str(width)))

figure(2),
plot(cte.fh, abs(TL.helmholtz1)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")

figure(3),
plot(cte.fh, abs(Z_HR)), xlabel("f [Hz]"), ylabel("TR - Helmholtz resonator [dB]")
    