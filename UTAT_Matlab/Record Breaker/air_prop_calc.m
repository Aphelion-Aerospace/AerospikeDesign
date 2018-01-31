%This function calculates air density and temperature at a specified altitude.
%Author: Abdullah Azzam
%Date: 10/31/2016

function [air_den,air_temp] = air_prop_calc(h)
% Input is altitude in metres
% Output is air density kg/m3, temp in K

%Parameters to be used

P0=101325; %P0: standard level atmospheric pressure, 101325 Pa
T0 = 288.15;%T0: sea level standard temperature, 288.15 K
g=9.81; %g: gravitational acceleration, 9.81 m/s2
L=0.0065;%L: temperature lapse rate, 0.0065 K/m
R=8.31447;%R: universal gas constant, 8.31447 J/(mol.K)
M=0.0289644;%M: molar mass of dry air, 0.0289644 kg/mol

%Temperature approximation inside troposphere at altitude h

T = T0 - (L*h);
air_temp = T;

%Pressure determination at altitude h

P = P0*(1-((L*h)/(T0))^((g*M)/(R*L)));

%Air density in kg/m3
air_den = (P*M)/(R*T); 

end
