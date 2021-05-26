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
mic.io.A = read_table(readtable('mic1e.csv','NumHeaderLines',1));
mic.io.B = read_table(readtable('mic2e.csv','NumHeaderLines',1));
mic.io.C = read_table(readtable('mic3e.csv','NumHeaderLines',1));
mic.io.D = read_table(readtable('mic4e.csv','NumHeaderLines',1));

% Partitioning the expansion chamber
mic.prt.A = read_table(readtable('mic1a.csv','NumHeaderLines',1));
mic.prt.B = read_table(readtable('mic2a.csv','NumHeaderLines',1));
mic.prt.C = read_table(readtable('mic3a.csv','NumHeaderLines',1));
mic.prt.D = read_table(readtable('mic4a.csv','NumHeaderLines',1));

% Vibro acoustics
mic.vibro.A = read_table(readtable('mic1v.csv','NumHeaderLines',1));
mic.vibro.B = read_table(readtable('mic2v.csv','NumHeaderLines',1));
mic.vibro.C = read_table(readtable('mic3v.csv','NumHeaderLines',1));
mic.vibro.D = read_table(readtable('mic4v.csv','NumHeaderLines',1));

% Absorbing material
mic.ab20.A = read_table(readtable('mic1ab20.csv','NumHeaderLines',1));
mic.ab20.B = read_table(readtable('mic2ab20.csv','NumHeaderLines',1));
mic.ab20.C = read_table(readtable('mic3ab20.csv','NumHeaderLines',1));
mic.ab20.D = read_table(readtable('mic4ab20.csv','NumHeaderLines',1));

mic.ab50.A = read_table(readtable('mic1ab.csv','NumHeaderLines',1));
mic.ab50.B = read_table(readtable('mic2ab.csv','NumHeaderLines',1));
mic.ab50.C = read_table(readtable('mic3ab.csv','NumHeaderLines',1));
mic.ab50.D = read_table(readtable('mic4ab.csv','NumHeaderLines',1));

mic.ab80.A = read_table(readtable('mic1ab80.csv','NumHeaderLines',1));
mic.ab80.B = read_table(readtable('mic2ab80.csv','NumHeaderLines',1));
mic.ab80.C = read_table(readtable('mic3ab80.csv','NumHeaderLines',1));
mic.ab80.D = read_table(readtable('mic4ab80.csv','NumHeaderLines',1));

% Perforated tube
mic.perf.A = read_table(readtable('mic1p.csv','NumHeaderLines',1));
mic.perf.B = read_table(readtable('mic2p.csv','NumHeaderLines',1));
mic.perf.C = read_table(readtable('mic3p.csv','NumHeaderLines',1));
mic.perf.D = read_table(readtable('mic4p.csv','NumHeaderLines',1));

mic.perf_small.A = read_table(readtable('mic1c.csv','NumHeaderLines',1));
mic.perf_small.B = read_table(readtable('mic2c.csv','NumHeaderLines',1));
mic.perf_small.C = read_table(readtable('mic3c.csv','NumHeaderLines',1));
mic.perf_small.D = read_table(readtable('mic4c.csv','NumHeaderLines',1));

mic.perf_large.A = read_table(readtable('mic1c1.csv','NumHeaderLines',1));
mic.perf_large.B = read_table(readtable('mic2c1.csv','NumHeaderLines',1));
mic.perf_large.C = read_table(readtable('mic3c1.csv','NumHeaderLines',1));
mic.perf_large.D = read_table(readtable('mic4c1.csv','NumHeaderLines',1));

% Perforated tube + absorbing material
mic.ab_tot.A = read_table(readtable('mic11.csv','NumHeaderLines',1));
mic.ab_tot.B = read_table(readtable('mic22.csv','NumHeaderLines',1));
mic.ab_tot.C = read_table(readtable('mic33.csv','NumHeaderLines',1));
mic.ab_tot.D = read_table(readtable('mic44.csv','NumHeaderLines',1));

% Helmholtz
mic.helm.A = read_table(readtable('mic1pa.csv','NumHeaderLines',1));
mic.helm.B = read_table(readtable('mic2pa.csv','NumHeaderLines',1));
mic.helm.C = read_table(readtable('mic3pa.csv','NumHeaderLines',1));
mic.helm.D = read_table(readtable('mic4pa.csv','NumHeaderLines',1));
%% Power 
power.normal.real = importdata('power_normal.csv').data(:,2);
power.normal.imag = importdata('power_normal.csv').data(:,3);

power.prt.real = importdata('vibro_divided.csv').data(:,2);
power.prt.imag = importdata('vibro_divided.csv').data(:,3);

power.expCh.real = importdata('vibro_sim-solution.csv').data(:,2);
power.expCh.imag = importdata('vibro_sim-solution.csv').data(:,3);

power.steel05.real = importdata('vibro0_5steel.csv').data(:,2);
power.steel05.imag = importdata('vibro0_5steel.csv').data(:,3);
power.steel15.real = importdata('vibro1_5steel.csv').data(:,2);
power.steel15.imag = importdata('vibro1_5steel.csv').data(:,3);
power.alum15.real = importdata('vibro1_5alum.csv').data(:,2);
power.alum15.imag = importdata('vibro1_5alum.csv').data(:,3);
power.tit05.real = importdata('vibro_tit05.csv').data(:,2);
power.tit05.imag = importdata('vibro_tit05.csv').data(:,3);

end

