%% To do list:
% add initial condition for launch angle
% add step incrementing angle of rocket, with no cross wind it is just
% going to be the launch angle

D = ___; %diameter of rocket
dtheta = 0; %change in the angle of rocket with respect to vertical, this is zero since we are not accoutning for crosswind yet
theta = zeros(1,_____); %initializing theta

if i==1
    theta(i) = a_launch;
end

[rho,T] = air_prop_calc(s(2,i)); %get air properties at given altitude
V = norm(v(:,i)); %get previous velocity magnitude
c = 331.3*sqrt(1+(T/273.15)); %calculate speed of sound at given temp
M = V/c; %calculate Mach number
C_D = coefficient_of_drag(GEOM, s(2,i), M, engine_mode); %get coefficient of drag
A = pi*D^2*0.25; %calculate cross sectional area of rocket, where D is the diameter
F_D = 0.5*rho*(v(:,i))^2*C_D*A; %calculate drag force vector
thrust_vec=[F(i)*sin(theta(i));F(i)*cos(theta(i))]; %get the thrust vector by projecting it along x and y given theta
a = ((thrust_vec - F_D)/m)-[0;0;g]; %get the acceleration vector
v(:,i+1)=v(:,i)+a*dt; %get velocity component at next time step
s(:,i+1)=s(:,i)+v(:,i)*dt+0.5*a*dt^2; %get position at next time step
theta(i+1) = theta(i) + dtheta; %changing angle of rocket at next time step