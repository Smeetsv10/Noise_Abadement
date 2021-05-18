function table = read_table(T)

table.freq = T.x_Frequency_Hz_;
table.p = T.Pressure_N_mm_2_MPa___Magnitude;
table.phase = T.Angle_degrees__Phase*(pi/180); 
end

