clear
close all

%% TMM
cte = set_cte();
mic = read_mics();
S_s = (pi/4)*cte.D^2;

vel.helm1.A = read_vel(readtable('vel2h.csv','NumHeaderLines',1));
vel.helm1.D = read_vel(readtable('vel3h.csv','NumHeaderLines',1));
vel.helm2.A = read_vel(readtable('vel2hv.csv','NumHeaderLines',1));
vel.helm2.D = read_vel(readtable('vel3hv.csv','NumHeaderLines',1));
vel.helm3.A = read_vel(readtable('vel2hv1.csv','NumHeaderLines',1));
vel.helm3.D = read_vel(readtable('vel3hv1.csv','NumHeaderLines',1));

for i = 1:length(cte.fh)   
    f = cte.fh(i);
    
    p1(i) = mic.helm1.A.p(i)*10^6;
    p2(i) = mic.helm1.D.p(i)*10^6;
    v1(i) = vel.helm1.A.v(i);
    v2(i) = -vel.helm1.D.v(i);
    R = [p2(i), p1(i); -v2(i), -v1(i)];
    P = [p1(i), p2(i); v1(i), v2(i)];
    TM = R/P;
    TM_11_mag(i) = (TM(1,1));
    TM_12_mag(i) = (TM(1,2));
    TM_21_mag(i) = (TM(2,1));
    TM_22_mag(i) = (TM(2,2));
    
    TL.HR(i) = 20*log10(abs((TM(1,1)+ (S_s/cte.c)*TM(1,2) + (cte.c/S_s)*TM(2,1) + TM(2,2))/2));
 
    
%     p_i = mic.helm3.A.p(i)*10^6;
%     v_xi = vel.helm3.A.v(i);
%     p_o =  mic.helm3.D.p(i)*10^6;
%     v_xo = vel.helm3.D.v(i);
%     
%     syms T11 T12 T21 T22
%     eq1 = T11*p_i+T12*v_xi == p_o;
%     eq2 = T11*p_o+T12*v_xo == p_i;
%     eq3 = T21*p_i+T22*v_xi == v_xo;
%     eq4 = T21*p_o+T22*v_xo == v_xi;
%     
%     [A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [T11, T12, T21, T22]);
%     T = linsolve(A,B);
% 
%     TM(i,:,:) = [T(1), T(2); T(3), T(4)];
%     TL.HR(i) = 20*log10( abs(T(1) + (S_s/cte.c)*T(2) + (cte.c/S_s)*T(3) + T(4)) / 2 );
end

%% Plotting
% figure(1)
% plot(cte.fh, TM_22_mag), xlabel('f [Hz]'), ylabel('T(1,1)')

figure(2),
plot(cte.fh, TL.HR), xlabel('f [Hz]'), ylabel('TL [dB]')

%% Function
function table = read_vel(T)
    cte = set_cte();
    S_s = (pi/4)*cte.D^2;
    
    table.freq = T.x_Frequency_Hz_;
    table.v = cte.rho_air*S_s*T.Velocity_mm_sec__Magnitude/10^-3;
    table.phase = T.Angle_degrees__Phase*(pi/180); 
end


