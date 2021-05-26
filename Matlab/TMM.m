clear
close all

%% TMM
mic = read_mics();
vel = readtable('air_cavity_sim2.csv','NumHeaderLines',1);
f = (20:2.5:2000);
for i = 1:length(f)
    if isnan(vel.Var1(i))
        vel(i,:) = [];
    end
end
disp(vel)

for i = 1:length(mic.expCh.B.p)
p_i = mic.expCh.B.p(i);
v_xi = 10;
p_o =  mic.expCh.C.p(i);
v_xo = 1.2;

syms T11 T12 T21 T22
eq1 = T11*p_i+T12*v_xi == p_o;
eq2 = T11*p_o+T12*v_xo == p_i;
eq3 = T21*p_i+T22*v_xi == v_xo;
eq4 = T21*p_o+T22*v_xo == v_xi;

[A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [T11, T12, T21, T22]);
T = linsolve(A,B);

TM = [T(1), T(2); T(3), T(4)];

end