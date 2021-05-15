clear

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
    close
    % Reflective muffler
    figure(4),
    subplot(3,1,1), hold on
    plot(cte.f, abs(TL.helmholtz1)), xlabel("f"), ylabel("TR - Helmholtz resonator [dB]")
    plot(cte.f, abs(TL.helmholtz2)), xlabel("f"), ylabel("TR - Helmholtz resonator [dB]")

    subplot(3,1,2)
    plot(cte.f, abs(TL.lambda4)), xlabel("f"), ylabel("TR - lambda/4 [dB]")
    subplot(3,1,3)
    plot(cte.f, abs(TL.expansion)), xlabel("f"), ylabel("TR - expansion chamber [dB]")
    
%     figure(5),
%     % Total TL:
%     plot(cte.f, abs(TL.expansion)+abs(TL.helmholtz)), xlabel("f"), ylabel("TR - expansion chamber")
    
    figure(6),
    plot(cte.f, abs(TL.expansion_NX), cte.f, abs(TL.expansion))
end

%save 'version1_2.mat'
