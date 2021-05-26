clear
close all

%% TMM
cte = set_cte();
mic = read_mics();

for i = 1:length(cte.fh)   
    f = cte.fh(i);

    
    velx1 = table2array(NXreader_PCH('h.pch').velocity_x.values);
    velx2 = NXreader_PCH('hv.pch').velocity_x;
    
    

%     for j = 1:length(mic.expCh.B.p)
%         p_i = mic.expCh.B.p(j);
%         v_xi = 10;
%         p_o =  mic.expCh.C.p(j);
%         v_xo = 1.2;
% 
%         syms T11 T12 T21 T22
%         eq1 = T11*p_i+T12*v_xi == p_o;
%         eq2 = T11*p_o+T12*v_xo == p_i;
%         eq3 = T21*p_i+T22*v_xi == v_xo;
%         eq4 = T21*p_o+T22*v_xo == v_xi;
% 
%         [A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [T11, T12, T21, T22]);
%         T = linsolve(A,B);
% 
%         TM = [T(1), T(2); T(3), T(4)];
%     end
end

%% Plotting
plot(cte.fh, velx1(1))


