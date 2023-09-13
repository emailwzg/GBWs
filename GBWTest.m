%{
Generalized Beta wavelets
Zhiguo Wang
Xi'an Jiaotong University
2023-09-113
cite: Wang, Zhiguo, Bing Zhang, Jinghuai Gao, Qingzhen Wang, and Qing Huo Liu, 
The wavelet transform with generalized Beta wavelets for seismic time-frequency analysis, 
Geophysics, 2017, 82(4), O47-O56
%}


clc
clear all
close all





n = 1024;
fs = 1000;
dt = 1/fs;
t = (-n/2+1:n/2)/fs;
t0 = 0.03;


a1 = 81;
b1 = 81;
a2 =81;
b2 = 3;

fm1 = 80;
fm2 = 50;

f11 = -0.5*(1-2*(pi*fm1*t).^2) .* exp(-(pi*fm1*t).^2);
f12 = -0.5*(1-2*(pi*fm2*t).^2) .* exp(-(pi*fm2*t).^2);
f1 = f11+f12;
f21 = 0.5*(1-2*(pi*fm1*(t-t0)).^2) .* exp(-(pi*fm1*(t-t0)).^2);
f22 = 0.5*(1-2*(pi*fm2*(t-t0)).^2) .* exp(-(pi*fm2*(t-t0)).^2);
f2 = f21+f22;
s = f1+f2;

fmax = 500;
fmin = 0.5 ;
      

[wave_gmw1,f0_gbw1] = gbwswavespdecomfun(s,a1,b1,fmin,fmax,dt);

[wave_gmw2,f0_gbw2] = gbwswavespdecomfun(s,a2,b2,fmin,fmax,dt);




figure(1)
subplot(411)
plot(t,s);
xlim([-0.1 0.1])
title('(a) Seismic wavelets')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(412)
pcolor(t,f0_gbw1,abs(wave_gmw1))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-0.1 0.1])
ylim([0 100])
title('(b) GBW, a=81,b=81')
shading interp
set(gcf,'color','w')

subplot(413)
pcolor(t,f0_gbw2,abs(wave_gmw2))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-0.1 0.1])
ylim([0 100])
title('(c) GBW, a=81,b=3')
shading interp
set(gcf,'color','w')









