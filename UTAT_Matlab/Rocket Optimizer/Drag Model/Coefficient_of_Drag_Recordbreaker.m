function [ Cd ] = Coefficient_of_Drag_Recordbreaker( Drag_Model, h, M, Rocket_Length, engine_mode)
%Coefficient_of_Drag 
% Description: This program interpolates the current coefficient of drag of a rocket based off it's altitude and mach number.
%  
% Written by Andreas Marquis and the propulsion gang - January 2017
%
%------------------

original_length = 2.911;
plus50_length = 1.5*original_length;

if h > 19812
    h = 19812;   %temporary fix
end

if (strcmpi(engine_mode, 'on'))
    
Cd_sea_level_original_length = interp1(Drag_Model(:,1),Drag_Model(:,3),M);
Cd_65k_original_length = interp1(Drag_Model(:,1),Drag_Model(:,5),M);

Cd_sea_level_plus50_length = interp1(Drag_Model(:,1),Drag_Model(:,7),M);
Cd_65k_plus50_length = interp1(Drag_Model(:,1),Drag_Model(:,9),M);

Cd_rocket_length_sea_level = interp1([original_length;plus50_length],[Cd_sea_level_original_length;Cd_sea_level_plus50_length],Rocket_Length);
Cd_rocket_length_65k = interp1([original_length;plus50_length],[Cd_65k_original_length;Cd_65k_plus50_length],Rocket_Length);

Cd = interp1([0;19812],[Cd_rocket_length_sea_level;Cd_rocket_length_65k],h);

elseif (strcmpi(engine_mode, 'off'))

Cd_sea_level_original_length = interp1(Drag_Model(:,1),Drag_Model(:,2),M);
Cd_65k_original_length = interp1(Drag_Model(:,1),Drag_Model(:,4),M);

Cd_sea_level_plus50_length = interp1(Drag_Model(:,1),Drag_Model(:,6),M);
Cd_65k_plus50_length = interp1(Drag_Model(:,1),Drag_Model(:,8),M);

Cd_rocket_length_sea_level = interp1([original_length;plus50_length],[Cd_sea_level_original_length;Cd_sea_level_plus50_length],Rocket_Length);
Cd_rocket_length_65k = interp1([original_length;plus50_length],[Cd_65k_original_length;Cd_65k_plus50_length],Rocket_Length);

Cd = interp1([0;19812],[Cd_rocket_length_sea_level;Cd_rocket_length_65k],h);
    
end

