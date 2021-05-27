clear
close all

%% Initialize parameters
cte = set_cte();
figures = [false, false, false, true, false]; % Expansion chamber, Vibro Acoustics, Optimisation, Absorbing, HR


%% Calculations
[IL,TL,ILh,TLh, power]= calculations_muffler(cte);

%% Post Processing
% Expansion chamber
if figures(1)
    
    figure(1),
    plot(cte.f, abs(TL.expansion), cte.f, abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('TL-Calculated', 'TL-Simulated')
    
    figure(2),
    plot(cte.f, IL.expansion, cte.f, IL.expansion_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('IL-Calculated', 'IL-Simulated')  
end

% Optimisation
if figures(2)
    
    figure(3),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.inlet_outlet_NX), cte.f, abs(TL.partitioning_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning')
    
    figure(4),
    plot(cte.f, IL.expansion_NX, cte.f, IL.inlet_outlet_NX, cte.f, IL.partitioning_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning') 
end

% Vibro acoustics
if figures(3)
    
    figure(5),
    plot(cte.f, power.expCh.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    
    figure(6),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.vibro_NX), cte.f, abs(TL.vibro5_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Vibro Acoustics steel 1,5mm','Vibro Acoustics steel 0,5mm')

    figure(7),
    plot(cte.f, IL.expansion_NX, cte.f, IL.vibro_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Vibro Acoustics') 

    figure(8),
    plot(cte.f, power.steel05.dB, cte.f, power.steel15.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Steel - 0.5mm', 'Steel - 1.5mm')

    figure(9),
    plot(cte.f, power.steel05.dB, cte.f, power.tit05.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Steel - 0.5mm', 'Titanium - 0.5mm')
    
    figure(10),
    plot(cte.f, power.expCh.dB, cte.f, power.prt.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Expansion chamber', 'Optimisations')
end

% Absorbing
if figures(4)
    
    % Absorbing material
    
    
    figure(11),
    plot(cte.f, abs(TL.ab20_NX), cte.f, abs(TL.ab50_NX), cte.f, abs(TL.ab80_NX),cte.f,abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL - Absorbing [dB]")
    legend('20mm', '50mm','80mm','Expansion chamber')
    figure(12),
    plot(cte.f, IL.ab20_NX, cte.f, IL.ab50_NX, cte.f, IL.ab80_NX, cte.f, IL.expansion_NX ), xlabel("f [Hz]"), ylabel("IL - Absorbing [dB]")
    legend('20mm', '50mm','80mm', 'Expansion chamber')
    
    % Perforated tube
    figure(),
    plot(cte.f, abs(TL.perf_NX), cte.f, abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL - Perforated [dB]")
    legend('with perforated', 'without perforated')
    
    figure(13),
    plot(cte.f, abs(TL.perf_small_NX), cte.f, abs(TL.perf_large_NX),cte.f,abs(TL.inlet_outlet_NX)), xlabel("f [Hz]"), ylabel("TL - Perforated [dB]")
    legend('small holes', 'large holes','inlet/outlet')
    figure(14),
    plot(cte.f, IL.perf_small_NX, cte.f, IL.perf_large_NX, cte.f, IL.inlet_outlet_NX), xlabel("f [Hz]"), ylabel("IL - Perforated [dB]")
    legend('small holes', 'large holes','inlet/outlet')
    
    % Final absorbing
    figure(15),
    plot(cte.f, abs(TL.ab_tot) ,cte.f,abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Absorbing', 'Expansion chamber')
    figure(16),
    plot(cte.f, IL.ab_tot, cte.f, IL.expansion_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Absorbing', 'Expansion chamber')
    
    % Comparison
    figure(17),
    plot(cte.f, abs(TL.ab_tot), cte.f, abs(TL.partitioning_NX), cte.f, abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Absorbing', 'Reflective','Expansion chamber')
    figure(18),
    plot(cte.f, IL.ab_tot, cte.f, IL.partitioning_NX, cte.f, IL.expansion_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Absorbing', 'Reflective','Expansion chamber')
end
    
if figures(5)

    figure(18),
    plot(cte.fh, abs(TLh.helm)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('TL-Calculated', 'TL-Simulated')
    

end
    


muffler_design(cte)
tilefigs

%save 'version1_2.mat'
