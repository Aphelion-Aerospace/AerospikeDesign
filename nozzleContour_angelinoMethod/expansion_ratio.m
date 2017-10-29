function epsilon = expansion_ratio(M,gamma)
	% calculates expansion ratio for gas expanding from choked conditions

	epsilon = 1./M.*(2./(gamma+1).*(1 + (gamma - 1)./2.*M.^2)).^((gamma+1)./(2.*(gamma-1)));

end
