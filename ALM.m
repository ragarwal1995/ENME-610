function lagrange_multiplier
clear all; close all


function gamma = func(x)
x1 = x(1); %Cathode thickness
x2 = x(2); %Anode thickness
x3 = x(3); %Cathode porosity
x4 = x(4); %Anode porosity
x5 = x(5); %Separator layer thickness
lambda = x(6);


%r_L = x_ccn + x2 + x_s + x1 + x_ccp
r_L = x2 + x5 + x1

%S = r_cell / r_L
S = 0.01 / r_L

%V_ccn = h_cell * (pi) * x_ccn * S * (r_L * (S-1) + x_ccn)
%V_ccn = 0.07 * (pi) * 10e-6 * S * (r_L * (S-1) + 10e-6)

%V_n = h_cell * (pi) * x2 * S * (r_L * (S-1) + x2 + 2*x_ccn)
V_n = 0.07 * (pi) * x2 * S * (r_L * (S-1) + x2)

%V_s = h_cell * (pi) * x_s * S * (r_L * (S-1) + x_s + 2*x_ccn + 2*x2)
V_s = 0.07 * (pi) * x5 * S * (r_L * (S-1) + x5 + 2*x2)

%V_p = h_cell * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x_ccn + 2*x2 + 2*x_s)
V_p = 0.07 * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x2 + 2*x5)

%V_ccp = h_cell * (pi) * x_ccp * S * (r_L * (S-1) + x_ccp + 2*x_ccn + 2*x2 + 2*x_s + 2*x1)
%V_ccp = 0.07 * (pi) * 15e-6 * S * (r_L * (S-1) + 15e-6 + 2*10e-6 + 2*x2 + 2*52e-6 + 2*x1)

%M_cell = M_ccn + M_n + M_s + M_p + M_ccp + M_e
%M_cell = (rho_ccn * V_ccn) + (rho_n * V_n * (1-x4)) + (rho_s * V_s * (1-epsilon_s)) + (rho_p * V_p * (1-x3)) + (rho_ccp * V_ccp) + rho_e * ((epsilon_s * V_s) + (x4 * V_n)+ (x3 * V_p))
M_cell =  (2270 * V_n * (1-x4)) + (900 * V_s * (1-0.46)) + (4140 * V_p * (1-x3)) + 1210 * ((0.46 * V_s) + (x4 * V_n)+ (x3 * V_p))

Q = 22055 * V_n * (1 - x4) * 96487
gamma  = (Q / M_cell) + lambda*(((22055*V_n*(1-x4))/(V_p*(1-x3)))-22860);
end

function dLambda = dfunc(x)
    dLambda = nan(size(x));

       h = 1e-3; % this is the step size used in the finite difference.
       for i=1:numel(x)
           dx=zeros(size(x));
           dx(i) = h;
           dLambda(i) = (func(x+dx)-func(x-dx))/(2*h);
       end
end

X1 = fsolve(@dfunc,[1 1 1 1 1 0],optimset('display','off'))
% let's check the value of the function at this solution
fval1 = func(X1)
sprintf('%f',fval1)
end


