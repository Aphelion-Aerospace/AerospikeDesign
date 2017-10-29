function [t_ratio,p_ratio,rho_ratio] = isen_stag_ratios(M,gamma)
%% takes in vector of Mach numbers, M, and ratio of specific heats, gamma
% returns the isentropic stagnation ratios (ie. T = T_0.*t_ratio)
isen_ratio = (1 + (gamma-1)./2.*M.^2);

t_ratio = 1./isen_ratio;

p_ratio = 1./(isen_ratio.^(gamma./(gamma-1)));

rho_ratio = 1./(isen_ratio.^(1./(gamma-1)));

end
