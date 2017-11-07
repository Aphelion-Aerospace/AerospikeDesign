function nu = prandtl_meyer(M,gamma)
% computes the prandtl meyer expansion angle for a given mach number, M, and ratio of specific heats, gamma

nu = sqrt((gamma+1)./(gamma-1)).*atan(sqrt((gamma-1)./(gamma+1).*(M.^2-1))) - atan(sqrt(M.^2-1));

end