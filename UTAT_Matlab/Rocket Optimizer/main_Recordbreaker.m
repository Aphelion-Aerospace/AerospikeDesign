%Main script for simulating full rocket performance. Uses engine simulation
%and 2D kinematics simulation to determine altitude of rocket. Employs
%Euler's Method to step variables through time. Uses pre-combustion values 
%in kinematics model to provide conservative estimate. Ignores ignition 
%transient and throat erosion. Plots monitored variables against time.
%Original scripts: Thomas S.H. Leung
%Revisions: Mitchell L. Passarelli
%Version Date: Saturday, March 11, 2017

%%%%%GENRAL NOTES:
%%%%%   In Thomas' file "Project14.m" I use full line comments to track
%%%%%   lines I have either accounted for or do not need in this file.

%%%%%NEW MEMBER TO-DO:
%%%%%   1.) Implement vapour burn.
%%%%%   2.) Implement cross-wind wind effects.

%%%%%SHORT-TERM TO-DO:
%%%%%   1.) See if I can use analytical expression for port diameter to
%%%%%   find burn time (by solving d_i - d(t) = 0) to get length of
%%%%%   matrices of monitored variables. Also see if I can couple it with
%%%%%   ox_tank model and pressure calculator to find the minumum of the
%%%%%   three burnout conditions.
%%%%%   2.) Test code by compaing with Thomas'. Make sure to check that the
%%%%%   results for the ox tank model are the same (esp. because of the use
%%%%%   of the n_l and n_g equations).

%%%%%LONG-TERM TO-DO (i.e., after finished):
%%%%%   1.) Look for a way to vectorize the loop (start with engine calcs).
%%%%%       Looks like only the ox tank model and fuel regression need to
%%%%%       be looped.
%%%%%   2.) Implement a check for choked flow and execute different
%%%%%       functions for choked & unchoked flow.

clear; clc; close all;

%% Import Inputs From Excel & Initialize Variables %%
%importing values into imported
% imported = import_params('Engine Preliminary Design Calculations.xlsx', 'Import Data - Recordbreaker', 1, 26);

%General parameters:
its = 2000;          %Est'd # iterations.
dt = 0.01;           %Time step size (s).
theta = 0;        %Launch angle, from vertical (deg).
rail_len = 15;              %Launch rail length (m).
comb_eff = 0.9;    %Combustive efficiency.
nozz_eff = 1;    %Nozzle efficiency.
P_a = 67461.44663;          %Ambient pressure at launch altitude (Pa).
T_a = 290.35;          %Ambient temperature at launch altitude (K).
m = 45;            %Rocket dry mass (kg).
rho_fuel = 950;     %Mass density of fuel (kg/m^3).
rocket_length = 6;   %Total length of rocket (m).

%CEA interpolation data:
num_rows = 356;     %# rows to import in CEA data file (= last row # - 1).
data = zeros(num_rows, 4);  %Cols. = OF, T, MW, gamma respectively.
[~, data(:, 1), data(:, 2), data(:, 3), ~, data(:, 4)] = ...
    import_cea_data('CEA_data.csv', 2, num_rows + 1);   %Get CEA data pts.
CEA = zeros(size(unique(data(:, 1)), 1), 4);    %Actual interp. pts matrix.
CEA(:, 1) = unique(data(:, 1));                 %Get unique OF values.

%Loop to clean up CEA data (average for each OF ratio):
for i = 1:size(CEA, 1)
    CEA(i, :) = mean(data(data(:, 1) == CEA(i, 1), :));
end
        %NOTE: Don't need this loop if you use interp2 with OF and Pcc.

%Initialize variables and constants:
i = 1;              %Loop index reset.
g0 = 9.81;          %Acceleration due to gravity at sea-level (m/s^2).
R_u = 8314;         %Universal gas constant (J/kmol*K).
MW_ox = 44.013;     %Molecular weight/mass of N2O (kg/kmol).
v_D = 0;            %Drag [velocity] loss (m/s).
f_D = 0.03;
L_fs = 0.00889;
d_fs = 0.0254;
A_inj = 9.60607*10.^-5;
Cd = 0.6;
s = zeros(2, its);  %Position vector (m).
v = s;              %Velocity vector (m/s).
a = s;              %Acceleration vector (m/s^2).
F = zeros(its, 1);  %Thrust (N).
Isp = F;            %Specific impulse (s).
OF = F;             %OF ratio.
Pcc = F;            %Combustion chamber pressure (Pa).
Pcc(1) = P_a;
Tcc = F;            %Combustion chamber temperature (K).
m_dot_ox = F;       %Ox mass flow rate (kg/s).
m_dot_f = F;        %Fuel mass flow rate (kg/s).
P_ox = F;           %Ox tank pressure array (Pa).
C_D = F;            %Drag coefficient array
F_D = F;            %Drag force array (N).
m_f = F;            %Fuel mass array (kg).
m_f(1) = 5.082836133;       %Initial fuel mass (kg).
diam = 0.1397;        %Diameter of rocket (m).
dtheta = 0;     %change in the angle of rocket with respect to vertical 
                %this is zero since we are not accoutning for crosswind yet    

%Ox tank parameters:
m_ox = F;
m_ox(1) = 22.8727626;      %Loaded liquid ox mass in tank (kg).
T_ox = T_a;                 %Initial ox tank temperature (K).
V_tank = 0.0;      %Internal volume of ox tank (m^3).
m_tank = 18.5578;      %Dry mass of ox tank (kg).
P_ox(1) = 4826330.105;      %Initial ox tank pressure (Pa).
Q = [2.781; 0.27244; 309.57; 0.2882];
        %Coefficients for molar specific vol. of liq. N2O (m^3/kmol).
Vhat_l = Q(2)^(1 + (1 - T_ox/Q(3))^Q(4))/Q(1);
        %Initial molar specific vol. of liq. N2O (m^3/kmol).
n(1) = m_ox(1)/MW_ox;       %# moles of liquid ox; n(2) is gaseous moles.
n_tot = n(1) + P_ox(1)*(V_tank - Vhat_l*n(1))/(R_u*T_ox);
n(2) = P_ox(1)*(V_tank - n_tot*Vhat_l)/(R_u*T_ox - P_ox(1)*Vhat_l);
        %NOTE: Assuming P_sat = P_ox(1).

%Fuel core parameters:
a0 = 0.000155;          %Regression rate coefficient (m/s^2).
n_reg = 0.5;    %FLUX EXP???   %Regression rate exponent.
num_ports = 1%%???  %# combustion ports in fuel grain.
L = 0.3048;         %% comb. chamberfull length??  %Fuel core length (m).
w = 0.0249;      %?????     %Initial fuel web thickness (m).
r_fo = 0.068/2          %Inner combustion chamber radius (m).
A_cc = pi*r_fo^2;           %Combustion chamber cross-section area (m^2).

%Nozzle parameters:
A_t = 0.0041;     %???PLug????    %[Initial] throat area (m^2).
A_e = 0.0048;         %Exit area (m^2).
e = @(Ma_e, k) (2/(k+1)*(1 + (k-1)/2*Ma_e^2))^((k+1)/(2*k-2))/Ma_e;
        %Function handle for nozzle area expansion ratio.
% e = @(Ma_cc, k) A_cc/A_e*Ma_cc/(1 + (k-1)/2*Ma_cc^2)^((k+1)/(2*k-2));
%         %Function handle for chamber-exit area-Mach relation.
nps = 3;
addpath('Drag Model')
Drag_Model = Recordbreaker_AeroLab_Drag_Model(nps);

%% Determine whether flow sys or injector plate primarily constricts flow:
%%%%%NOTE: Don't worry about this one yet. Realistically, the injector
%%%%%plate should always constrict the flow.


%% Simulation Time Loop (ignition to apogee) %%
%Loop until rocket reaches apogee (ignore condition on first iteration):
while v(2, i) > 0 || i == 1
    engine_mode = 'off';
    
    %If burn not finished, continue engine simulation:
    if m_f(i) > 0 && Pcc(i) <= P_ox(i) && m_ox(i) > 0
        burnout_alt = s(2, i);      %Store current alt. as burnout.
        engine_mode = 'on';
        
        %Simulate blow-down in ox tank:
        if i == 1   %Uses initial guesses for initial Pcc and m_dot_ox.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n, T_ox, V_tank, m_tank, 2.8, f_D, L_fs, d_fs, ...
                A_inj, Cd, 3.44e6);
        else        %Uses prev. iteration data.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n, T_ox, V_tank, m_tank, m_dot_ox(i - 1), f_D, ...
                L_fs, d_fs, A_inj, Cd, Pcc(i - 1));
        end
        
        %Determine mass flow rates and fuel regression:
        r_port = r_fo - w;                  %Get fuel port radius.
        G_ox = m_dot_ox(i)/(pi*r_port^2);   %Calc. ox mass flux.
        reg_rate = a0*G_ox^n_reg;           %Calc. fuel regression rate.
        SA_port = 2*pi*r_port*L;            %Calc. fuel core surface area.
        m_dot_f(i) = reg_rate*rho_fuel*SA_port*num_ports;
        m_dot_prop = m_dot_ox(i) + m_dot_f(i);  %Total mass flow rate.
        
        %Determine propellant thermochemical properties:
        OF(i) = m_dot_ox(i)/m_dot_f(i);
        Tcc(i) = interp1(CEA(:,1), CEA(:,2), OF(i));
        gamma = interp1(CEA(:,1), CEA(:,4), OF(i));
        R = R_u/interp1(CEA(:,1), CEA(:,3), OF(i));
        char_vel = comb_eff*sqrt(gamma*R*Tcc(i))/gamma * ...
            ((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
        
        %Eval. engine performance including nozzle:
        if i > 10   %Calc. combustion chamber P (avg to damp out osc.).
            Pcc(i) = (m_dot_prop*char_vel/A_t + Pcc(i - 1))/2;
        else
            Pcc(i) = m_dot_prop*char_vel/A_t;
        end
        Ma_e = fzero(@(Ma) A_e/A_t - e(Ma, gamma), 2.75);   %Exit Mach.
        P_e = Pcc(i)/(1 + (gamma-1)/2*Ma_e^2)^(gamma/(gamma-1));
                %Flow is always choked for Pcc >= 200 psi.
        
%         %Eval. engine performance when flow not choked:
%         rho = Pcc(i)/R/Tcc(i);              %Calc. chamber gas density.
%         v_cc = m_dot_prop/rho/A_cc;         %Get flow velocity.
%         Ma_cc = v_cc/sqrt(gamma*R*Tcc(i));  %Combustion chamber Mach #.
%         Ma_e = fzero(@(Ma_e) Ma_e/(1 + (gamma-1)/2*Ma_e)^...
%             ((gamma+1)/(2*gamma-2)) - e(Ma_cc, gamma), 0.5);    %Exit Mach.
        %%%%%Don't worry about above until sim is finished: can implement a
        %%%%%check for choked condition: look at AER310 midterm Q3 to see
        %%%%%how to check for choking.
        
        v_e = sqrt(2*gamma*R*Tcc(i)/(gamma - 1) * ...
            (1 - (P_e/Pcc(i))^((gamma - 1)/gamma)));    %Get exit velocity.
        F(i) = nozz_eff*(m_dot_prop*v_e + (P_e - P_a)*A_e);     %Thrust.
        Isp(i) = F(i)/m_dot_prop/g0;        %Calc. specific impulse.
        
        %Update engine parameters:
        %Make sure to store all parameters in (i+1)th entry, but pass only
        %the ith entries to the kinematics model.
        w = w - reg_rate*dt;                %Regress fuel core.
        n = n + dn*dt;                      %Consume moles of ox.
        m_f(i + 1) = rho_fuel*L*pi*w*(2*r_fo - w);  %Update fuel mass.
        m_ox(i + 1) = n(1)*MW_ox;           %Update ox mass.
%         m_f(i + 1) = m_f(i) - m_dot_f(i)*dt;        %Consume fuel mass.        
%         m_ox(i + 1) = m_ox(i) + dn(1)*MW_ox*dt;     %Consume ox mass.
        T_ox = T_ox + dT_ox*dt;             %Step ox tank temp.
        %%%%Include updates to ox tank parameters here too.%%%%
        
    %Change here for tracking
    else
        Pcc(i) = P_a;
        P_ox(i) = P_a;
        m_f(i + 1) = m_f(i);
        if m_f(i) > 0
            m_ox(i + 1) = 0;
        else
            m_ox(i+1) = m_ox(i);
        end
        
        %%%%%blow down model goes here!!!%%%%%
        
    end

    %Move rocket forward in time using kinematics model:
    %%%%Call to kinematics model. Note that this needs to output the
    %%%%ambient pressure at the current altitude into P_a.
    %%%%Pull apart Thomas' KinematicsP14 function into atomic blocks and
    %%%%code my script here to use those blocks and pass them into the
    %%%%kinematics model.
    %%%%PSEUDOCODE:
    %%%%    1.) Calc. ambient air properties at current altitude.
    %%%%    2.) Calc. drag based on current conditions.
    %%%%    3.) Pass drag coefficient, etc. into dynamics/kinematics model.
    %%%% ------------------------------------------------------------------
    [rho, T] = air_prop_calc(s(2,i));  %get air properties at given alt.
    V = norm(v(:,i));                %get previous velocity magnitude
    c = sqrt(1.4*287*T);             %calc. speed of sound at given temp
    M = V/c;                         %calculate Mach number
    
    if M > 3    M = 3;  end
    
    C_D(i) = Coefficient_of_Drag_Recordbreaker(Drag_Model, s(2,i), M, ...
        rocket_length, engine_mode);          %get coefficient of drag
              
    S = pi/4*diam^2;	%Calc. cross-sectional area of rocket.
    F_D(i) = rho/2*(v(1, i)^2 + v(2, i)^2)*S*C_D(i);    %Drag mag..
    
    %Check if rocket is still on launch rail and change angle accordingly:
    if s(2, i) > rail_len
            theta = atan2d(v(1, i), v(2, i));
    end
    
    a(:, i) = (F(i) - F_D(i))*[sind(theta); cosd(theta)] ...
        /(m + m_f(i) + m_ox(i)) - [0; g0];          %Calc. accel.
    v(:, i+1) = v(:, i) + a(:, i)*dt;               %Calc. new velocity.
    s(:, i+1) = s(:, i) + v(:, i)*dt + a(:, i)/2*dt^2;	%Get new position.
    
    v_D = v_D + F_D(i)*cosd(theta)/(m + m_f(i) + m_ox(i))*dt;   %Accum. drag loss.
	%%%% ------------------------------------------------------------------
    
    i = i + 1;      %Increment counter index.
end

%Trim monitored variables:
m_dot_ox = m_dot_ox(1:i);       %ox mass flow rate
F = F(1:i);                     %Thrust
Pcc = Pcc(1:i);                 %Combustion Chamber Pressure
P_ox = P_ox(1:i);               %Ox tank pressure
Tcc = Tcc(1:i);                 %Combustion chamber temperature
OF = OF(1:i);                   %Ox/Fuel mass ratio
C_D = C_D(1:i);                 %Drag Coef
F_D = F_D(1:i);                 %Drag force
v = v(:, 1:i);                  %Velocity componentsf
s = s(:, 1:i);                  %Position components
m_f = m_f(1:i);                 %Fuel mass
m_ox = m_ox(1:i);               %ox mass
t = [0:dt:(i - 1)*dt];          %Time vector.

%% Plot Outputs/Monitored Variables %%
%close all
figure(1)
subplot(2,1,1);
plot(t, s(1, :), t, s(2, :),...
    'linewidth',2);
legend('x direction','y direction', 'Location', 'NW')
grid on
ylabel('Position')

subplot(2,1,2);
plot(t, v(1, :), t, v(2, :),...
    'linewidth',2);
legend('x direction','y direction')
grid on
xlabel('time [s]')
ylabel('Velocity')

figure(2)
subplot(3,1,1);
plot(t,m_ox,t,m_f,...
    'linewidth',2);
legend('m_{ox} [kg]','m_{fuel} [kg]')
grid on
xlabel('time [s]')
ylabel('Mass [kg]')

subplot(3,1,2);
plot(t,m_dot_ox,...
    'linewidth',2);
grid on
xlabel('time [s]')
ylabel('Rate of m_{ox}')

subplot(3,1,3);
plot(t,OF,...
    'linewidth',2);
grid on
xlabel('time [s]')
ylabel('O/F ratio')

figure(3)
subplot(2,1,1);
plot(t,Pcc,t,P_ox,...
    'linewidth',2);
legend('Combustion Chamber Pressure','Ox Tank Pressure')
grid on
xlabel('time [s]')
ylabel('Pressure [Pa]')

subplot(2,1,2);
plot(t,Tcc,...
    'linewidth',2);
grid on
xlabel('time [s]')
ylabel('Flame Temperature')

figure(4)
subplot(3,1,1);
plot(t,F,...
    'linewidth',2);
grid on
ylabel('Thrust [N]')

subplot(3,1,2);
plot(t,F_D,...
    'linewidth',2);
grid on
xlabel('time [s]')
ylabel('Drag Force [N]')

subplot(3,1,3);
plot(t,C_D,...
    'linewidth',2);
grid on
xlabel('time [s]')
ylabel('Drag Coefficient')

Export_array = [v_D;  burnout_alt; s(2,end)]
%fprintf('Recordbreaker_output.txt',Export_array);

%Export velocity and position data (for flutter)
for i = 1:numel(v(1,:))
    vel_norm(i) = norm(v(:,i));
end

Export_flutter_data = [s(2, :)',  vel_norm'];
xlswrite('Rocket_Positionv_v_Speed.xlsx',Export_flutter_data, 'Export Data');

