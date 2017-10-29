function [x,r,M] = angelino_nozzle_contour(r_e,r_b,M_e, gamma,points)
% calculates contour of an aerospike engine using Angelino's method

% given the exit radius, r_e, end of plug radius, r_b, desired exit mach
% number, M_e, ratio of specific heats, gamma, and the number of contour
% points to be approximated

% returns length (x) and radius (r) pairs that define the axisymmetric
% nozzle contour. Also returns the mach number vector used to define
% contour. Note, this method assumes the nozzle lip is located at (0,0)

% estimating contour
A_e = pi.*(r_e.^2-r_b.^2);
A_t = A_e./expansion_ratio(M_e,gamma);

M = linspace(1,M_e,points);

A = A_t.*expansion_ratio(M,gamma);

alpha = prandtl_meyer(M_e,gamma) - prandtl_meyer(M,gamma) + mach_angle(M);
l = (r_e-sqrt(r_e.^2 - (A.*M.*sin(alpha)./pi)))./sin(alpha);

x = l.*cos(alpha);
r = l.*sin(alpha);

end