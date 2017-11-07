function mach_zero_func = mach_expansion_ratio(M,gamma,epsilon)

mach_zero_func = 1./M.*(2./(gamma+1).*(1 + (gamma - 1)./2.*M.^2)).^((gamma+1)./(2.*(gamma-1))) - epsilon;

end