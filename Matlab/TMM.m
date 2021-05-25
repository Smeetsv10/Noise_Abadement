clear
close all

%% TMM
p_i = 10;
v_xi = 10;
p_o = 10;
v_xo = 1.2;

syms T11 T12 T21 T22
eq1 = T11*p_i+T12*v_xi == p_o;
eq2 = T11*p_o+T12*v_xo == p_i;
eq3 = T21*p_i+T22*v_xi == v_xo;
eq4 = T21*p_o+T22*v_xo == v_xi;

[A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [T11, T12, T21, T22]);
T = linsolve(A,B);

TM = [T(1), T(2); T(3), T(4)]