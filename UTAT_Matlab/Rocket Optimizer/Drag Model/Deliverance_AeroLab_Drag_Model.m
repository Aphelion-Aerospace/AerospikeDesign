% Deliverance_AeroLab_Drag_Model.m
%
% author: Andreas Marquis
% Date: Winter 2017
%
% Description: This program estimates the drag on a rocket. Based off the
% the drag estimation software: AeroLab 1.3.2.
% 
%           By Andreas Marquis
%
%
%Aerolab is based on the following litterature:
% 
% USAF Stability and control DATCOM
%         Flight Control Division
%         Air Force Flight Dynamics Laboratory
%         Wright-Patterson Air Force Base, Ohio
%         1978
% 
% Boundary Layer Theory
%         Dr. H Schlichting
%         McGraw-Hill
% 
% Mechanics of Fluids
%         W. J. Duncan, A .S Thom, A. D. Young
%         Edward Arnold 
% 
% Missile Configuration Design
%         S. S. Chin
%         McGraw-Hill 1961
% 
% NACA Tecnical Note 1032
% NACA Technical Note 4201
% NACA Tecnical Report 1386
% NACA Tecnical Report 1153
% NACA Tecnical Report 1353
% NACA Tecnical Report 1374
% NACA Tecnical Report 1253
% NACA Tecnical Report 1227
% NACA Tecnical Report 1238
% NACA Tecnical Report 1307
% NACA Research Memorandum RM A53A30
% 
% Nasa Nike Tomahawk Handbook
%         Thiokol Chemical Corporation
%         Astro-Met Division
%         P.O.Box 1497 Ogden, Utah
%         1966
% 
% Stabilisering af raketter ved Hjælp af Barrowmans metode
%         J. Franck
%         Dansk Amatør Raket Klub
% 
% Theory of Wing Sections
%         Ira h. Abbot, A. E. Von Doenhoff
%         Dover Publications, Inc 1958
% 
% MIL-HDBK-762(MI): DESIGN OF AERODYNAMICALLY STABILIZED FREE ROCKETS
% 
% All values of the atmospheric parameters are from the US STANDARD ATMOSPHERE 
% SUPPLEMENTS, 1966 (60deg N, July)
%
% History
%-------------------
% 01-15-17 Created by AM

function [  Drag_Model ] = Deliverance_AeroLab_Drag_Model()

Drag_Model = NaN(101,5);

dat = xlsread('Deliverance II Drag Curve','sea level');
Drag_Model(1:end,1) = dat(1:end,1); %importing mach numbers 
Drag_Model(1:end,2) = dat(1:end,2); %importing unpowered drag @ ground level
Drag_Model(1:end,3) = dat(1:end,3); %importing powered drag @ ground level

dat = xlsread('Deliverance II Drag Curve','30,000 feet');
Drag_Model(1:end,4) = dat(1:end,2); %importing unpowered drag @ 30,000 feet
Drag_Model(1:end,5) = dat(1:end,3); %importing powered drag @ 30,000 feet

end
