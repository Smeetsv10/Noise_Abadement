function [IL, TL] = calc_IL_TL(mic, cte, i)
f = cte.f(i);
k = 2*pi*f/cte.c;

A1 = (mic.A.p(i)*exp(1i*mic.A.phase(i))*exp(1i*k*cte.s1)-mic.B.p(i)*exp(1i*mic.B.phase(i)))*exp(-1i*k*(cte.m1+cte.s1))/(exp(1i*k*cte.s1)-exp(-1i*k*cte.s1)); % amp inlet
A3 = (mic.C.p(i)*exp(1i*mic.C.phase(i))*exp(1i*k*cte.s2)-mic.D.p(i)*exp(1i*mic.D.phase(i)))*exp(-1i*k*(cte.m2-cte.d))/(exp(1i*k*cte.s2)-exp(-1i*k*cte.s2)); % amp outlet
%A3 = mic.C.p(i);

IL = 20*log10(abs(cte.A3_wo/A3));
TL = 10*log10(abs(A1/A3)^2); %10*log(abs((mic_A.p(i)/mic_C.p(i))^2))
    
end

