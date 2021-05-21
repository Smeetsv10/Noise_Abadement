function [IL,TL] = absorption_muffler(IL, TL, power, cte)
%% Absorption type muffler
    %% Import Data
    [mic, power] = read_mics();

    %% Parameters
    p_i = 10; % incident pressure
    p_t = 1; % transmitted pressure
    
    %% TPP muffler
    TL.TPP = 20*log((S_1/S_2)^0.5*abs(p_i/p_t));
    
    [IL.absorption_NX(i),TL.absorption_NX(i)] = calc_IL_TL(mic.abs, cte, i);
end

