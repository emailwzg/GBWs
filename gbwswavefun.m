%{
Generalized Beta wavelets
Zhiguo Wang
Xi'an Jiaotong University
2023-09-13
cite: Wang, Zhiguo, Bing Zhang, Jinghuai Gao, Qingzhen Wang, and Qing Huo Liu, 
The wavelet transform with generalized Beta wavelets for seismic time-frequency analysis, 
Geophysics, 2017, 82(4), O47-O56
%}
function W=gbwswavefun(k,a,b,scale);
n = length(k);
f= 2*atan(scale.*k)/pi;
% alpha=1/beta(2*a+1,2*b+1)
alpha=1;
expnt = (1-f).^b.*(k > 0.);  
norm1=  sqrt(scale*k(2))*sqrt(n)*alpha*(f).^a;
daughter = norm1.*expnt;
daughter1 = daughter.*(k > 0.);
W =daughter1/norm(daughter1,2);
% W=daughter;


