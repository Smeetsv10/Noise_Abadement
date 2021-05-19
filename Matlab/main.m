clear
close

%% Initialize parameters
cte = set_cte();
figures = true;

%% Calculations
% Reflective Type muffler
[IL,TL]= reflective_muffler(cte);

% Absorption type muffler

% Meta muffler


%% Post Processing
if figures
    
    % Reflective muffler
    figure(4),
    subplot(4,1,1), hold on
    plot(cte.f, abs(TL.helmholtz1)), xlabel("f"), ylabel("TR - Helmholtz resonator [dB]")
    plot(cte.f, abs(TL.helmholtz2)), xlabel("f"), ylabel("TR - Helmholtz resonator [dB]")
    subplot(4,1,2)
    plot(cte.f, abs(TL.lambda4)), xlabel("f"), ylabel("TR - lambda/4 [dB]")
    subplot(4,1,3)
    plot(cte.f, abs(TL.expansion)), xlabel("f"), ylabel("TR - expansion chamber [dB]")
    subplot(4,1,4), hold on
    plot(cte.f, abs(TL.helmholtz1)), xlabel("f"),
    plot(cte.f, abs(TL.helmholtz2)), xlabel("f"),
    plot(cte.f, abs(TL.expansion)), xlabel("f"),
    plot(cte.f, abs(TL.lambda4)), xlabel("f"), ylabel("TR [dB]")
    legend("Helmholtz resonator 1","Helmholtz resonator 2", "expansion chamber", "lambda/4 ")
        
    figure(5),
    TL.total = TL.expansion + TL.lambda4;
    plot(cte.f, abs(TL.total)), xlabel("f"), ylabel("total TR [dB]")
    
    figure(6),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.expansion)), xlabel("f"), ylabel("TL [dB]")
    legend('TL-Calculated', 'TL-Simulated')
    
    figure(7),
    plot(cte.f, IL.expansion,cte.f, IL.expansion_NX), xlabel("f"), ylabel("IL [dB]")
    legend('IL-Calculated', 'IL-Simulated')
        
end

%save 'version1_2.mat'
