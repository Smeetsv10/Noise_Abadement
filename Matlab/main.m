clear
close all
%

%% Initialize parameters
cte = set_cte();
figures = [false, true, false, false]; % Total, Expansion chamber, Optimisation, Vibro Acoustics

%% Calculations
[IL,TL, power]= calculations_muffler(cte);

%% Post Processing
% Total reflective
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

% Expansion chamber
if figures(2)
    
    figure(3),
    plot(cte.f, abs(TL.expansion), cte.f, abs(TL.expansion_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('TL-Calculated', 'TL-Simulated')
    
    figure(4),
    plot(cte.f, IL.expansion, cte.f, IL.expansion_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('IL-Calculated', 'IL-Simulated')  
end

% Optimisation
if figures(3)
    
    figure(5),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.inlet_outlet_NX), cte.f, abs(TL.partitioning_NX)), xlabel("f [Hz]"), ylabel("TL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning')
    
    figure(6),
    plot(cte.f, IL.expansion_NX, cte.f, IL.inlet_outlet_NX, cte.f, IL.partitioning_NX), xlabel("f [Hz]"), ylabel("IL [dB]")
    legend('Expansion chamber', 'Inlet/outlet', 'Partitioning') 
end

% Vibro acoustics
if figures(4)
    
    figure(7),
    plot(cte.f, power.expansion.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    
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
    plot(cte.f, power.steel15.dB, cte.f, power.alum15.dB), xlabel("f [Hz]"), ylabel("Radiated power [dB]")
    legend('Steel - 1.5mm', 'Aluminium - 1.5mm')

    figure(12),
    plot(cte.f, power.steel15.dB, cte.f, power.alum15.dB), xlabel("thickness [mm]"), ylabel("Price [eur]")
    legend('Steel - 1.5mm', 'Aluminium - 1.5mm')


end

muffler_design(cte)
tilefigs

%save 'version1_2.mat'
