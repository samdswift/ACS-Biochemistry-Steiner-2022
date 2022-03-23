clear
clc

%[R-C-T7] - scan
%[P-N-T7] - 10-fold R-C-T7. The quick graph is below from Matlab:

%code (using equation 7 with P = 10xR): 
R = logspace(0,4,200); 
Kd2 = 2; 
Kd1 = 1000;
EC50 = 0.5.*R + Kd2./19 + Kd1./(9.*R./Kd2 - 1);
loglog(R,EC50)

writematrix(transpose(R),'5fx.csv')
writematrix(transpose(EC50),'5fy.csv')