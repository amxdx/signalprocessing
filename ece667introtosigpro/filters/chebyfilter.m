%Design a Butterworth LowPass filter for the given specifications 
%Given Ordash, Ou , K1, K2, Design LPF HPF Cheby
%Ordash = 430 rad/sec; Ou = 300 rad/sec; K1 = -3db; K2 = -25db;  
clear all; clc;
Ordash = 430; 
Ou = 300; 
Or = Ordash/Ou; 
K1 = -3; 
K2 = -25; 
A = 10^(-K2/20); 
E = sqrt(10^(-K1/10) - 1 ) ; 
g = sqrt( (A^2 - 1)/ E^2); 
n = ceil(log10(g + sqrt(g^2 -1))/log10( Or + sqrt(Or^2 -1)));
for i = 1:3
[z,p,k] = buttap(n-1+i);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
sys = tf(num,den);            % Create transfer function 
[numt, dent] = lp2lp(num,den,Ou);
[h,w] = freqs(numt,dent);
mag=abs(h);
phase=angle(h);
subplot(2,1,1);
plot(w,mag,'DisplayName',['Order '  num2str(n-1+i)''])
legend('show')
axis([0 1000 0 1]);
title('Magnitude Response')
grid on
hold on
subplot(2,1,2);
plot(w,phase,'DisplayName',['Order '  num2str(n-1+i)''])
grid on
axis([0 1000 1.5*min(phase) 1.5*max(phase)]);
title('Phase Response')
hold on
end
suptitle(['Order ' num2str(n)' ', ' num2str(n+1) 'and ' num2str(n+2)' 'Low-Pass Butterworth Filter'])
figure; 
for i = 1:3
[z,p,k] = buttap(n-1+i);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
sys = tf(num,den);            % Create transfer function 
[numt, dent] = lp2lp(num,den,Ou);
[h,w] = freqs(numt,dent);
mag=abs(h);
phase=angle(h);
subplot(2,1,1);
plot(w,mag,'DisplayName',['Order '  num2str(n-1+i)''])
legend('show')
axis([0 1000 0 1]);
title('Magnitude Response')
grid on
hold on
subplot(2,1,2);
plot(w,phase,'DisplayName',['Order '  num2str(n-1+i)''])
grid on
axis([0 1000 1.5*min(phase) 1.5*max(phase)]);
title('Phase Response')
hold on
end
suptitle(['Order ' num2str(n)' ', ' num2str(n+1) 'and ' num2str(n+2)' 'High-Pass Butterworth Filter'])