rho_fuel = 950 
L = 0.3048 #fuel core length CHECK
r_fo = 0.068/2 # Inner combustion chamber radius
a0 = 0.000155 # Regression rate coeff (m/s**2)
n_reg = 0.5 # Regression rate exponent, FLUX EXP???
MW_ox = 44.013 # Molecular weight/mass of N2O (kg/kmol)
m_ox = #liquid ox mass in tank initial

n = m_ox/MW_ox # moels of liquid ox

r_port = r_fo - w #fuel port radius; r_fo = inner combustion chamber radius

G_ox = m_dot_ox[i]/(np.pi*r_port**2) #ox mass flux

reg_rate = a0*G_ox**n_reg # n_reg = reg. rate exp.; a0 = reg. rate coeff.

w = w - reg_rate*dt # w = initial fuel web thickness

n = n + dn*dt # consume moles of ox

m_f[i+1] = rho_fuel*L*np.pi*w*(2*r_fo-w) # mass fuel