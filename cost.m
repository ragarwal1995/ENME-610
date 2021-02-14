function f0_hat = cost(x,Rp)
x1 = x(1); %Cathode thickness
x2 = x(2); %Anode thickness
x3 = x(3); %Cathode porosity
x4 = x(4); %Anode porosity
x5 = x(5); %Separator layer thickness
%r_L = x_ccn + x2 + x_s + x1 + x_ccp
r_L = x2 + x5 + x1;

%S = r_cell / r_L
S = floor(0.01 / r_L);

%V_ccn = h_cell * (pi) * x_ccn * S * (r_L * (S-1) + x_ccn)
%V_ccn = 0.07 * (pi) * 10e-6 * S * (r_L * (S-1) + 10e-6)

%V_n = h_cell * (pi) * x2 * S * (r_L * (S-1) + x2 + 2*x_ccn)
V_n = 0.07 * (pi) * x2 * S * (r_L * (S-1) + x2);

%V_s = h_cell * (pi) * x_s * S * (r_L * (S-1) + x_s + 2*x_ccn + 2*x2)
V_s = 0.07 * (pi) * x5 * S * (r_L * (S-1) + x5 + 2*x2);

%V_p = h_cell * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x_ccn + 2*x2 + 2*x_s)
V_p = 0.07 * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x2 + 2*x5);

M_cell =  (2270 * V_n * (1-x4)) + (900 * V_s * (1-0.46)) + (4140 * V_p * (1-x3)) + 1210 * ((0.46 * V_s) + (x4 * V_n)+ (x3 * V_p));

Q = 22055 * V_n * (1 - x4) * 96487;

C = ((22055*V_n*(1-x4))/((V_p*(1-x3))))-22860;

ubx1 = 250e-6;
ubx2 = 250e-6;
ubx3 = 0.6;
ubx4 = 0.6;
ubx5 = 100e-6;

lbx1 = 40e-6;
lbx2 = 40e-6;
lbx3 = 0.2;
lbx4 = 0.2;
lbx5 = 10e-6;

C1 = x1-ubx1;
C2 = x2-ubx2;
C3 = x3-ubx3;
C4 = x4-ubx4;
C5 = x5-ubx5;
C6 = -x1+lbx1;
C7 = -x2+lbx2;
C8 = -x3+lbx3;
C9 = -x4+lbx4;
C10 = -x5+lbx5;

f = -(Q / M_cell);

f0_hat = f + Rp*((max(0,C)^2)+(max(0,C1)^2)+(max(0,C2)^2)+(max(0,C3)^2)+(max(0,C4)^2)+(max(0,C5)^2)+(max(0,C6)^2)+(max(0,C7)^2)+(max(0,C8)^2)+(max(0,C9)^2)+(max(0,C10)^2));