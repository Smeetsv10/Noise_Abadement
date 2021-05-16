function muffler_design(cte)

% coordinates:
inlet = [0 137];
ExpCh = inlet(2); % change if extended inlet
outlet = [ExpCh+cte.L*10^3 400];

mic_A = ExpCh-(cte.m1+cte.s1)*10^3;
mic_B = ExpCh-cte.m1*10^3;
mic_C = ExpCh+(cte.L+cte.m2)*10^3;
mic_D = ExpCh+(cte.L+cte.m2+cte.s2)*10^3;
HR_1 = 36;
HR_2 = 88;


figure(1),title("Muffler design"), hold on
plot(inlet,1000*[-cte.D/2 -cte.D/2],'b-'),
plot(inlet,1000*[cte.D/2 cte.D/2],'b-'),
plot([inlet(1) inlet(1)],1000*[-cte.D/2 cte.D/2],'b-'),

plot([ExpCh ExpCh],1000*[cte.D/2 cte.D*5/2],'b-'),
plot([ExpCh ExpCh],1000*[-cte.D/2 -cte.D*5/2],'b-'),
plot([ExpCh ExpCh+cte.L*10^3],1000*[-cte.D*5/2 -cte.D*5/2],'b-'),
plot([ExpCh ExpCh+cte.L*10^3],1000*[cte.D*5/2 cte.D*5/2],'b-'),
plot([ExpCh+cte.L*10^3 ExpCh+cte.L*10^3],1000*[cte.D/2 cte.D*5/2],'b-'),
plot([ExpCh+cte.L*10^3 ExpCh+cte.L*10^3],1000*[-cte.D/2 -cte.D*5/2],'b-'),

plot(outlet,1000*[-cte.D/2 -cte.D/2],'b-'),
plot(outlet,1000*[cte.D/2 cte.D/2],'b-'),
plot([outlet(2) outlet(2)],1000*[-cte.D/2 cte.D/2],'b-'),

plot([mic_A mic_B mic_C mic_D],[0 0 0 0],'ro'),

axis([-5 405 -205 205])

end

