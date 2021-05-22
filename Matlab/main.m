clear
close all

%% Initialize parameters
cte = set_cte();
figures = [false, false, false, true]; % Expansion chamber, Vibro Acoustics,  Optimisation, Absorbing


%% Calculations
[IL,TL, power]= calculations_muffler(cte);

%% Post Processing
% Expansion chamber
if figures(1)
    
    figure(3),
    plot(cte.f, abs(TL.expansion), cte.f, abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('TL-Calculated', 'TL-Simulated')
    
    figure(4),
    plot(cte.f, IL.expansion, cte.f, IL.expansion_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('IL-Calculated', 'IL-Simulated')  
end

% Optimisation
if figures(2)
    
    figure(5),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.inlet_outlet_NX), cte.f, abs(TL.partitioning_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning')
    
    figure(6),
    plot(cte.f, IL.expansion_NX, cte.f, IL.inlet_outlet_NX, cte.f, IL.partitioning_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning') 
end

% Vibro acoustics
if figures(3)
    
    figure(7),
    plot(cte.f, power.expCh.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    
    figure(8),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.vibro_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Vibro Acoustics')

    figure(9),
    plot(cte.f, IL.expansion_NX, cte.f, IL.vibro_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Vibro Acoustics') 

    figure(10),
    plot(cte.f, power.steel05.dB, cte.f, power.steel15.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Steel - 0.5mm', 'Steel - 1.5mm')

    figure(11),
    plot(cte.f, power.steel05.dB, cte.f, power.tit05.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Steel - 0.5mm', 'Titanium - 0.5mm')
    
    figure(12),
    plot(cte.f, power.expCh.dB, cte.f, power.prt.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Expansion chamber', 'Optimisations')
end

% Absorbing
if figures(4)
    
    figure(12),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.ab50_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Absorbing - 50mm')
    
    figure(13),
    plot(cte.f, IL.expansion_NX, cte.f, IL.ab50_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Absorbing -50mm') 
end

muffler_design(cte)
tilefigs

%save 'version1_2.mat'
