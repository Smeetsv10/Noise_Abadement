function [mic, power] = read_mics()
%% Microphones    
% Normal duct
mic.wo.A = read_table(readtable('mic1n','NumHeaderLines',1));
mic.wo.B = read_table(readtable('mic2n','NumHeaderLines',1));
mic.wo.C = read_table(readtable('mic3n','NumHeaderLines',1));
mic.wo.D = read_table(readtable('mic4n','NumHeaderLines',1));

% Expansion tube
mic.expCh.A = read_table(readtable('mic1.csv','NumHeaderLines',1));
mic.expCh.B = read_table(readtable('mic2.csv','NumHeaderLines',1));
mic.expCh.C = read_table(readtable('mic3.csv','NumHeaderLines',1));
mic.expCh.D = read_table(readtable('mic4.csv','NumHeaderLines',1));

% Inlet/outlet extension
mic.io.A = read_table(readtable('mic1a.csv','NumHeaderLines',1));
mic.io.B = read_table(readtable('mic2a.csv','NumHeaderLines',1));
mic.io.C = read_table(readtable('mic3a.csv','NumHeaderLines',1));
mic.io.D = read_table(readtable('mic4a.csv','NumHeaderLines',1));

% Partitioning the expansion chamber
mic.prt.A = read_table(readtable('mic1e1.csv','NumHeaderLines',1));
mic.prt.B = read_table(readtable('mic2e1.csv','NumHeaderLines',1));
mic.prt.C = read_table(readtable('mic3e1.csv','NumHeaderLines',1));
mic.prt.D = read_table(readtable('mic4e1.csv','NumHeaderLines',1));

% Vibro acoustics
mic.vibro.A = read_table(readtable('mic1v.csv','NumHeaderLines',1));
mic.vibro.B = read_table(readtable('mic2v.csv','NumHeaderLines',1));
mic.vibro.C = read_table(readtable('mic3v.csv','NumHeaderLines',1));
mic.vibro.D = read_table(readtable('mic4v.csv','NumHeaderLines',1));

%% Power 
power.normal.real = importdata('power_normal.csv').data(:,2);
power.normal.imag = importdata('power_normal.csv').data(:,3);

power.expansion.real = importdata('vibro_sim-solution.csv').data(:,2);
power.expansion.imag = importdata('vibro_sim-solution.csv').data(:,3);

power.steel.real(:,1) = importdata('vibro0_5steel.csv').data(:,2);
power.steel.imag(:,1) = importdata('vibro0_5steel.csv').data(:,3);
power.steel.real(:,2) = importdata('vibro1_5steel.csv').data(:,2);
power.steel.imag(:,2) = importdata('vibro1_5steel.csv').data(:,3);

power.alum.real = importdata('vibro1_5alum.csv').data(:,2);
power.alum.imag = importdata('vibro1_5alum.csv').data(:,3);

end

