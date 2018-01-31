function [ Cd ] = Coefficient_of_Drag( Drag_Model, h, M, rocket_length, engine_mode)
%Coefficient_of_Drag 
% Description: This program interpolates the current coefficient of drag of a rocket based off it's altitude and mach number.
%  
% Written by Andreas Marquis and the propulsion gang - January 2017
%
%------------------

if h > 9144
    h = 9144;   %temporary fix
end

if (strcmpi(engine_mode, 'on'))
    
Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,3),M);
Cd_30k = interp1(Drag_Model(:,1),Drag_Model(:,5),M);

Cd = interp1([0;9144],[Cd_sea_level;Cd_30k],h); % 9144 meters is 30,000 feet
    
elseif (strcmpi(engine_mode, 'off'))

Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,2),M);
Cd_30k = interp1(Drag_Model(:,1),Drag_Model(:,4),M);

Cd = interp1([0;9144],[Cd_sea_level;Cd_30k],h); % 9144 meters is 30,000 feet
    
end

