%alpha_range = [[0:0.1:0.5]1 3 5];
gamma = 2;
Rp = 1;
f_old = 0;
format long
x0 = [100e-6, 100e-6, 0.3, 0.3,55e-6]; %Initial point 
while Rp<1e9%||diff==0
    %alpha = alpha_range(counter);
    Rp = Rp * gamma;
    anonymousFunctionHandle = @(x) cost(x0,Rp);
    [xhatstar,f] = fminsearch(anonymousFunctionHandle,x0,optimset('Tolx',1e-10,'MaxFunEvals',10000,'MaxIter',10000))
    diff = f-f_old;
    f_old = f;
end