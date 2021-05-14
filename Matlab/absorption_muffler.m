function [IL,TL] = absorption_muffler(cte)
%% Absorption type muffler
    %% Import Data
    
    %% Parameters
    p_i = 10; % incident pressure
    p_t = 1; % transmitted pressure
    
    %% TPP muffler
    
    TL.TPP = 20*log((S_1/S_2)^0.5*abs(p_i/p_t));
    
end

