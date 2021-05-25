function [IL,TL,ILh,TLh, power] = calculations_muffler(cte)
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
    
    %% Normal
    [IL.normal_NX(i),TL.normal_NX(i)] = calc_IL_TL(mic.wo, cte, i);
    %% Inlet/outlet extension
    [IL.inlet_outlet_NX(i),TL.inlet_outlet_NX(i)] = calc_IL_TL(mic.io, cte, i);

    %% Partitioning expansion chamber
    [IL.partitioning_NX(i),TL.partitioning_NX(i)] = calc_IL_TL(mic.prt, cte, i);
    power.prt.abs(i) = sqrt(power.prt.real(i)^2+power.prt.imag(i)^2);
    power.prt.dB(i) = 10*log10(power.prt.abs(i) / cte.power_ref);
    
    %% Vibro acoustics
    power.expCh.abs(i) = sqrt(power.expCh.real(i)^2+power.expCh.imag(i)^2);
    power.normal.abs(i) = sqrt(power.normal.real(i)^2+power.normal.imag(i)^2);
    power.expCh.dB(i) = 10*log10(power.expCh.abs(i) / cte.power_ref);
 
    [IL.vibro_NX(i),TL.vibro_NX(i)] = calc_IL_TL(mic.vibro, cte, i);

    power.steel05.dB(i) = 10*log10( sqrt(power.steel05.real(i)^2+power.steel05.imag(i)^2) / cte.power_ref);
    power.steel15.dB(i) = 10*log10( sqrt(power.steel15.real(i)^2+power.steel15.imag(i)^2) / cte.power_ref);
    power.alum15.dB(i) = 10*log10( sqrt(power.alum15.real(i)^2+power.alum15.imag(i)^2) / cte.power_ref);
    power.tit05.dB(i) = 10*log10( sqrt(power.tit05.real(i)^2+power.tit05.imag(i)^2) / cte.power_ref);
    
    t(i) = (0.5 + 0.0050*i)*10^-3; %0.5..1.5
    Vol_mat(i) = 2*(pi/4)*(D_new^2-cte.D^2)*t(i) + 2*pi*D_new*(cte.L-2*t(i))*t(i);
    tot_steel_price(i) = cte.price_steel*cte.rho_steel*Vol_mat(i);
    tot_titanium_price(i) = cte.price_titanium*cte.rho_titanium*Vol_mat(i);
%% 
% Absorption type Mufflers   
% _________________________

    %% Absorbing material
    [IL.ab20_NX(i),TL.ab20_NX(i)] = calc_IL_TL(mic.ab20, cte, i);
    [IL.ab50_NX(i),TL.ab50_NX(i)] = calc_IL_TL(mic.ab50, cte, i);
    [IL.ab80_NX(i),TL.ab80_NX(i)] = calc_IL_TL(mic.ab80, cte, i);  
    
    %% Perforated ducts
    [IL.perf_NX(i),TL.perf_NX(i)] = calc_IL_TL(mic.perf, cte, i);
    [IL.perf_small_NX(i),TL.perf_small_NX(i)] = calc_IL_TL(mic.perf_small, cte, i);
    [IL.perf_large_NX(i),TL.perf_large_NX(i)] = calc_IL_TL(mic.perf_large, cte, i);

    %% Total
    [IL.ab_tot(i),TL.ab_tot(i)] = calc_IL_TL(mic.ab_tot, cte, i);
    

    
%% 
% Rest
% _________________________

    
end

for i = 1:length(cte.fh)
    %% Helmholtz
    
    [ILh.helm(i),TLh.helm(i)] = TL_helmholtz(mic.helm, cte, i);
end
%     figure(12),
%     plot(t*1000, tot_steel_price, t*1000, tot_titanium_price), xlabel("thickness [mm]"), ylabel("Price [eur]")
%     legend('Steel', 'Titanium')

end

