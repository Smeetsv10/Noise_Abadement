function [mic, vibro] = read_mics()
    
% Normal duct
mic.A_wo = read_table(readtable('mic1n','NumHeaderLines',1));
mic.B_wo = read_table(readtable('mic2n','NumHeaderLines',1));
mic.C_wo = read_table(readtable('mic3n','NumHeaderLines',1));
mic.D_wo = read_table(readtable('mic4n','NumHeaderLines',1));

% Expansion tube
mic.A = read_table(readtable('mic1.csv','NumHeaderLines',1));
mic.B = read_table(readtable('mic2.csv','NumHeaderLines',1));
mic.C = read_table(readtable('mic3.csv','NumHeaderLines',1));
mic.D = read_table(readtable('mic4.csv','NumHeaderLines',1));

% Inlet/outlet extension
mic.A_io = read_table(readtable('mic1a.csv','NumHeaderLines',1));
mic.B_io = read_table(readtable('mic2a.csv','NumHeaderLines',1));
mic.C_io = read_table(readtable('mic3a.csv','NumHeaderLines',1));
mic.D_io = read_table(readtable('mic4a.csv','NumHeaderLines',1));

% Partitioning the expansion chamber
mic.A_prt = read_table(readtable('mic1e1.csv','NumHeaderLines',1));
mic.B_prt = read_table(readtable('mic2e1.csv','NumHeaderLines',1));
mic.C_prt = read_table(readtable('mic3e1.csv','NumHeaderLines',1));
mic.D_prt = read_table(readtable('mic4e1.csv','NumHeaderLines',1));

% Vibro acoustics
vibro.data = importdata('vibro_sim-solution.csv').data;
end

