clear
clc

data = [
3000	0.382   0.106646873
3000	0.123   0.106646873
3000	0.282   0.106646873
1000	0.249   0.117575508
1000	0.105   0.117575508
1000	0.393   0.117575508
333     0.189   0.038938698
333     0.095   0.038938698
333     0.128   0.038938698
111     0.159   0.040745825
111     0.068   0.040745825
111     0.078   0.040745825
37      0.12    0.026545558
37      0.07    0.026545558
37      0.059   0.026545558
12.3	0.092   0.020066556
12.3	0.046   0.020066556
12.3	0.054   0.020066556
4.1     0.074   0.013123346
4.1     0.044   0.013123346
4.1     0.049   0.013123346]; 

Protein = data(:,1);
A = data(:,2);

fun = @(F) (F(1)*(Protein/(Protein + F(2)))-A);
[bfs] = lsqnonlin(fun,[20 15000]);
S = bfs(1)
K = bfs(2)
%%%%%

%S
answer = zeros(1000,2);
for i = 1:0.001:2
index = round(i*1000);
answer(index,1) = i;    
holder = zeros(length(data),1);    
    for j = 1:length(data)
    EST = i*data(j,1)/(data(j,1) + K);
    X2 = (((data(j,2)-0.044)-EST)./data(j,3)).^2;
    holder(j) = X2;
    end
answer(index,2) = sum(holder);
end

%K
answer = zeros(100000,2);
for i = 1:length(answer)
answer(i,1) = i;    
holder = zeros(length(data),1);    
    for j = 1:length(data)
    EST = S*data(j,1)/(data(j,1) + i);
    X2 = ((data(j,2)-0.044-EST)./data(j,3)).^2;
    holder(j) = X2;
    end
answer(i,2) = sum(holder);
end

%%%%%
semilogx(answer(:,1),answer(:,2))
