function muffler_design(cte)

% x-coordinates:
inlet = 0;
outlet = 400;
mic_A = 5;
mic_B = 15;
mic_C = 385;
mic_D = 395;
HR_1 = 36;
HR_2 = 88;
ExpCh = 137;

figure, hold on
plot([inlet outlet],1000*[-cte.D/2 -cte.D/2],'b-'),
plot([inlet HR_2-1000*cte.D*0.3],1000*[cte.D/2 cte.D/2],'b-'),

axis([-5 405 -205 205])

end

