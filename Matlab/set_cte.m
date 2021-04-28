function cte = set_cte()
% General parameters
cte.g = 9.81; % m/s^2
cte.c = 341; % m/s
cte.rho_air = 1.224; % density [kg/m^3]
cte.mu_air = 1.8*10^-5; % viscosity [kg/m*s]
cte.gamma = 1.4;

% Simcenter parameters [m, Hz]
cte.L = 0.300; % length
cte.H = 0.200; % height
cte.W = 0.200; % width
cte.D = 0.040; % diameter duct
cte.t = 0.001; %0.5..1.5
cte.f_step = 10;
cte.f = (20:cte.f_step:2000)';

% External parameters
cte.steel_price = 40; % eur/kg
cte.titanium_price = 80; % eur/kg


end

