clear all
clc
%% Define the parameters
global D Q1 Q2 Ct Kc K1 K2 T_inf h L1 L2 R1 R2
D = 20e-3; %Diameter of the cell in millimeter
Q1 = 1; %Heat generation in cathod W/m^3
Q2 = 3; %Heat generation in anode W/m^3
Ct = 20e-6; %Seperator layer thickness in micro meter
Kc = 20; %Seperator layer thermal conductivity
K1 = 10; %Cathode thermal conductivity
K2 = 5; %Anode thermal conductivity
L1 = 5.1e-3; %Cathode mean free path at 0% porosity
L2 = 3e-3; %Anode mean free path at 0% porosity
R1 = 1e-8; %Cathode particle size
R2 = 2e-9; %Anode particle size
T_inf = 20; %Coolant temperature in C
h = 10; %Coolant convection coeffiecent
%% Optimize
function Thermal
x0 = [100, 100, 0.3, 0.3]; %Initial point correct  
%sigma = 10^5; 
lb = [40, 40, 0.2, 0.2]; %Lower bound correct 
ub = [250, 250, 0.6, 0.6]; %Upper bound correct 
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%options = optimset('LargeScale', 'off');
[x, f] = fmincon(@Thermal_Objective, x0, [], [], [], [], lb, ub, @Thermal_Constraint, options) 
end

function f = Thermal_Objective(x)
x1 = x(1); %Cathode thickness
x2 = x(2); %Anode thickness
x3 = x(3); %Cathode porosity
x4 = x(4); %Anode porosity
S = D./(x1+x2-2*Ct); %Define S, the number of layers factor correct
T1 = (S*(Q1*x1+Q2*x2)/h)+T_inf; %Surface temperature
K1p = K1*(1-x3)./(1+(L1/R1)*x3.^(1/3));%Cathode thermal conductivity at x3 porosity
K2p = K2*(1-x4)./(1+(L2/R2)*x4.^(1/3));%Anode thermal conductivity at x4 porosity
Lambda = (x1+x2+Ct)./((x1./K1p)+(x2./K2p)+(Ct/Kc));%Equvelat thermal conductivity of the system
f = (S*(Q1+Q2)./(2*Lambda))*D^2+T1
%f = ((D/(x1+x2-2*Ct))*(Q1+Q2)/(2*((x1+x2+Ct)/((x1/(K1*(1-x3)/(1+(L1/R1)*x3^(1/3))))+(x2/(K2*(1-x4)/(1+(L2/R2)*x4^(1/3))))+(Ct/Kc)))))*(D^2)+(((D/(x1+x2-2*Ct))*(Q1*x1+Q2*x2)/h)+T_inf)
end

function [C,Ceq] = Thermal_Constraint(x)
x1 = x(1); %Cathode thickness
x2 = x(2); %Anode thickness
x3 = x(3); %Cathode porosity
x4 = x(4); %Anode porosity
C(1) = [];
C(2) = [];
Ceq = []; 
end