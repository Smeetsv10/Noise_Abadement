function cte = set_cte()
% General parameters
cte.g = 9.81; % m/s^2
cte.c = 341; % m/s
cte.rho_air = 1.224; % density [kg/m^3]
cte.mu_air = 1.8*10^-5; % viscosity [kg/m*s]
cte.gamma = 1.4;

% Simcenter parameters [m, Hz]
cte.L = 0.213; % length
cte.H = 0.200; % height
cte.W = 0.200; % width
cte.D = 0.040; % diameter duct
cte.t = 0.001; %0.5..1.5
cte.f_step = 10;
cte.f = (20:cte.f_step:2000)';
cte.m1 = 0.0785;
cte.s1 = 0.010;
cte.m2 = 0.0785;
cte.s2 = 0.010;
cte.d = 0.0935;

% External parameters
cte.price_steel = 2; % eur/kg
cte.rho_steel = 8050; %kg/m^3
cte.price_titanium = 45; % eur/kg
cte.rho_titanium = 4510; % eur/kg

end

