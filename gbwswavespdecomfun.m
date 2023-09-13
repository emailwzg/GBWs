%{
Generalized Beta wavelets
Zhiguo Wang
Xi'an Jiaotong University
2023-09-113
cite: Wang, Zhiguo, Bing Zhang, Jinghuai Gao, Qingzhen Wang, and Qing Huo Liu, 
The wavelet transform with generalized Beta wavelets for seismic time-frequency analysis, 
Geophysics, 2017, 82(4), O47-O56
%}

function  [wave,f0]=gbwswavespdecomfun(Y,a,b,fmin,fmax,dt);

n = length(Y);

 %c = tan(a*pi/2/(a+b));
 c = tan((2*a+1)*pi/2/(2*a+2*b+1));

s0 = c/(2*pi*fmax);
sJ = c/(2*pi*fmin);


dj = 1./68.;

J1=fix((log(sJ/s0)/log(2))/dj);
k = 1:fix(n/2);
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];
f = fft(Y);  
scale = s0*2.^((0:J1)*dj);
wave = zeros(J1+1,n); 
wave = wave + 1i*wave;
for a1 = 1:J1+1
	daughter=gbwswavefun(k,a,b,scale(a1));	
	wave(a1,:) = ifft(f.*daughter);  
end


f0 = c./(2*pi*scale);
