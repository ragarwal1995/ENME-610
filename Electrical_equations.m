clear all
clc
format long
%% Define the parameters
global D Q1 Q2 Kc K1 K2 T_inf h L1 L2 R1 R2 h_cell C_0_anode F rho_n rho_s rho_p rho_e rho_ccn rho_ccp x_ccn x_ccp  r_cell epsilon_s
%Energy density parameters
h_cell = 0.07; % height of cell in m
C_0_anode = 22055; %Concentraion of anode in mol/m^3
F = 96487; %Faraday's Constant in C/mol
rho_n = 2270; %Density of anode in kg/m^3
rho_s = 900; %Density of separator in kg/m^3
rho_p = 4140; %Density of cathode in kg/m^3
rho_e = 1210; %Density of electrolyte in kg/m^3
rho_ccn = 8700; % Density of negative current collector Cu foil kg/m^3
rho_ccp = 2700; % Density of positive current collector Al in foil kg/m^2
x_ccn = 0; % Thickness of negative current collector Cu foil in m (ASSUMED ZERO FOR SIMPLICITY)
x_ccp = 0; % Thick of positive current collector Al foil in m (ASSUMED ZERO FOR SIMPLICITY)
r_cell = 0.01; %m
epsilon_s = 0.46; %Porosity of Separator
%% Optimize

Energy
function Energy
x0 = [100e-6, 100e-6, 0.3, 0.3,55e-6]; %Initial point 
lb = [40e-6, 40e-6, 0.2, 0.2,10e-6]; %Lower bound 
ub = [250e-6, 250e-6, 0.6, 0.6,100e-6]; %Upper bound  

options = optimset('LargeScale', 'off');
%options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp');
[x, f] = fmincon(@Energy_Objective, x0, [], [], [], [], lb, ub, @nonlcon , options) 
end

function f = Energy_Objective(x)
    global D Q1 Q2 Kc K1 K2 T_inf h L1 L2 R1 R2 h_cell C_0_anode F rho_n rho_s rho_p rho_e rho_ccn rho_ccp x_ccn x_ccp  r_cell epsilon_s
    x1 = x(1); %Cathode thickness
    x2 = x(2); %Anode thickness
    x3 = x(3); %Cathode porosity
    x4 = x(4); %Anode porosity
    x5 = x(5); %Separator layer thickness
    r_L = x2 + x5 + x1;
    S = floor(r_cell / r_L)
    V_ccn = h_cell * (pi) * x_ccn * S * (r_L * (S-1) + x_ccn);
    V_n = h_cell * (pi) * x2 * S * (r_L * (S-1) + x2 + 2*x_ccn) %?
    V_s = h_cell * (pi) * x5 * S * (r_L * (S-1) + x5 + 2*x_ccn + 2*x2); %???
    V_p = h_cell * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x_ccn + 2*x2 + 2*x5); %?
    V_ccp = h_cell * (pi) * x_ccp * S * (r_L * (S-1) + x_ccp + 2*x_ccn + 2*x2 + 2*x5 + 2*x1);
    M_cell = (rho_ccn * V_ccn) + (rho_n * V_n * (1-x4)) + (rho_s * V_s * (1-epsilon_s)) + (rho_p * V_p * (1-x3)) + (rho_ccp * V_ccp) + rho_e * ((epsilon_s * V_s) + (x4 * V_n)+ (x3 * V_p))
    Q = C_0_anode * V_n * (1 - x4) * F
    f = -(Q / M_cell);
end

function [C,Ceq] = nonlcon(x)
    global D Q1 Q2 Kc K1 K2 T_inf h L1 L2 R1 R2 h_cell C_0_anode F rho_n rho_s rho_p rho_e rho_ccn rho_ccp x_ccn x_ccp  r_cell epsilon_s
    x1 = x(1); %Cathode thickness
    x2 = x(2); %Anode thickness
    x3 = x(3); %Cathode porosity
    x4 = x(4); %Anode porosity
    x5 = x(5); %Separator layer thickness
    r_L = x2 + x5 + x1;
    S = floor(r_cell / r_L);
    V_n = h_cell * (pi) * x2 * S * (r_L * (S-1) + x2 + 2*x_ccn) 
    V_p = h_cell * (pi) * x1 * S * (r_L * (S-1) + x1 + 2*x_ccn + 2*x2 + 2*x5); 
    C(1) = ((C_0_anode*V_n*(1-x4))/(V_p*(1-x3)))-22860; 
       Ceq = []; 
end