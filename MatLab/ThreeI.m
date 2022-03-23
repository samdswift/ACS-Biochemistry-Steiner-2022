clear
clc

KdD = 99999999999999999999999999999999999999999999999999999999999999; 
Kdoff = 10000; 
P = 1000;
answer = zeros(10000,1);
for R = 1:10000
    %%All values above are in nanoomolar
    p = [(2./KdD) (1 + 2.*Kdoff./KdD) (P + Kdoff - R) (-Kdoff.*R)];
    n = roots(p); Rmin = real(n(3));
    Cmin = Rmin.*P./(Kdoff + Rmin);
    nano = Cmin./1000
    answer(R) = nano;
    %answer[R,2] = Cmin
end
writematrix(answer,'Inf.csv')