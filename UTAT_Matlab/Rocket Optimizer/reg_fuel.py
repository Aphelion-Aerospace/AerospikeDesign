r_port = r_fo - w #fuel port radius; r_fo = inner combustion chamber radius

G_ox = m_dot_ox[i]/(np.pi*r_port**2) #ox mass flux

reg_rate = a0*G_ox**n_reg # n_reg = reg. rate exp.; a0 = reg. rate coeff.

w = w - reg_rate*dt # w = initial fuel web thickness

n = n + dn*dt # consume moles of ox

m_f[i+1] = rho_fuel*L*np.pi*w*(2*r_fo-w)