function muffler_design(cte)

% coordinates:
ExpCh = [-106.5 106.5];
inlet = [ExpCh(1)-50 ExpCh(1)];
outlet = [ExpCh(2) ExpCh(1)+50];

mic_A = ExpCh(1)-(cte.m1+cte.s1)*10^3;
mic_B = ExpCh(1)-cte.m1*10^3;
mic_C = ExpCh(2)+(cte.L+cte.m2)*10^3;
mic_D = ExpCh(2)+(cte.L+cte.m2+cte.s2)*10^3;
HR_1 = 36;
HR_2 = 88;


figure(1),title("Muffler design"), hold on
plot(inlet,1000*[-cte.D/2 -cte.D/2],'b-'),
plot(inlet,1000*[cte.D/2 cte.D/2],'b-'),
plot([inlet(1) inlet(1)],1000*[-cte.D/2 cte.D/2],'b-'),

plot([ExpCh(1) ExpCh(1)],1000*[cte.D/2 cte.D*5/2],'b-'),
plot([ExpCh(1) ExpCh(1)],1000*[-cte.D/2 -cte.D*5/2],'b-'),
plot(ExpCh,1000*[-cte.D*5/2 -cte.D*5/2],'b-'),
plot(ExpCh,1000*[cte.D*5/2 cte.D*5/2],'b-'),
plot([ExpCh(2) ExpCh(2)],1000*[cte.D/2 cte.D*5/2],'b-'),
plot([ExpCh(2) ExpCh(2)],1000*[-cte.D/2 -cte.D*5/2],'b-'),
% 
% plot(outlet,1000*[-cte.D/2 -cte.D/2],'b-'),
% plot(outlet,1000*[cte.D/2 cte.D/2],'b-'),
% plot([outlet(2) outlet(2)],1000*[-cte.D/2 cte.D/2],'b-'),

plot([mic_A mic_B mic_C mic_D],[0 0 0 0],'ro'),

axis([-200 200 -205 205])

end

