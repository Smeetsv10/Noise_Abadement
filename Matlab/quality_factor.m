function [Q,R_s] = quality_factor(cte, L, A, V, f, v)
 w = 2*pi*f;
 BL_thickness = sqrt(2*cte.mu_air/(cte.rho_air*w));
 if BL_thickness > cte.t/2
     h = BL_thickness;
 else
     h = cte.t/2;
 end
 eps = 0.5;
 M = v/cte.c;
 
 R_s = (cte.rho_air*cte.c/S_2)*((k*BL_thickness*(2*pi*D_new)*L/(2*S_2))*(1+(cte.gamma-1)*sqrt(5/(3*gamma))) ... % acoustic resistance tube
        + 0.288*k*BL_thickness*log(4*S_2/(pi*h^2)) + eps*S_2*k^2/(2*pi) + M);
   
 Q = (cte.rho_air*cte.c/R_s)*sqrt(L/(A*V)); 
end

