r_e = 0.15; % exit radius
r_b = 0.05; % plub nozzle base radius
M_e = 2.5; % desired exit mach number
T_0 = 350; % chamber temperature(300 K)
p_0 = 2e7; % ~200 bar
gamma = 1.4; % ratio of specific heats

[x,r,M] = angelino_nozzle_contour(r_e,r_b,M_e, gamma,100);
% note, angelino's method assumes nozzle lip is at (0,0)

plot(x,r,'k-',x,r_e.*ones(size(x)),'k--');
axis equal

legend('Contour','Axisymmetric Centre')

[t_ratio,p_ratio,rho_ratio] = isen_stag_ratios(M,gamma);
T = t_ratio.*T_0;
P = p_ratio.*p_0;