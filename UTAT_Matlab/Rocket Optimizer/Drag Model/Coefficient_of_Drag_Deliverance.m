function [C_D] = Coefficient_of_Drag_Deliverance(Drag_Model, h, M, burning)
%COEFFICIENT_OF_DRAG Interpolates to find current drag coefficient.
%   Returns drag coefficient of a rocket based on its current altitude and 
%   Mach # by interpolating the data specified by the input drag model.
%   INPUTS:
%       "Drag_Model" = ___________
%       "h" = scalar for current altitude [m].
%       "M" = scalar current Mach #.
%       "burning" = scalar logical where TRUE => engine is burning.
%   OUTPUT:
%       "C_D" = scalar drag coefficient.
% Written by Andreas Marquis and the propulsion gang - January 2017
% Revised: Friday, March 17, 2017 by Mitchell Passarelli
%------------------

if h > 9144
    h = 9144;   %temporary fix
end

if burning      %Reduce base drag while engine is burning.
    Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,3),M);
    Cd_30k = interp1(Drag_Model(:,1),Drag_Model(:,5),M);
    
    C_D = interp1([0; 9144],[Cd_sea_level; Cd_30k], h); %9144 m = 30,000 ft.
else            %Use normal base drag after burnout.
    Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,2),M);
    Cd_30k = interp1(Drag_Model(:,1),Drag_Model(:,4),M);
    
    C_D = interp1([0;9144],[Cd_sea_level;Cd_30k],h); %9144 m = 30,000 ft.
end

end
