clear

%% Initialize parameters
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
    subplot(3,1,1)
    plot(cte.f, abs(TL.helmholtz)), xlabel("f"), ylabel("TR - Helmholtz resonator")
    subplot(3,1,2)
    plot(cte.f, abs(TL.lambda4)), xlabel("f"), ylabel("TR - lambda/4")
    subplot(3,1,3)
    plot(cte.f, abs(TL.expansion)), xlabel("f"), ylabel("TR - expansion chamber")
    
%     figure(5),
%     % Total TL:
%     plot(cte.f, abs(TL.expansion)+abs(TL.helmholtz)), xlabel("f"), ylabel("TR - expansion chamber")
    
    figure(6),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.expansion))
end
