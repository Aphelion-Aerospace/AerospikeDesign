%% Design of nozzle + estimation of flow properties
r_e = 0.0762; % exit radius, check if this is diam. or radius
A_t = 0.004; % plug nozzle throat area
epsilon = 7; %expansion ratio
T_0 = 2527.77; % chamber temperature(300 K)
p_0 = 19.098; % ~200 bar
gamma = 1.83; % ratio of specific heats

[x,r,M] = angelino_nozzle_contour(A_t,epsilon, r_e,gamma,100);

lip_coord = [0;0]; % note, angelino's method assumes nozzle lip is at (0,0)



[t_ratio,p_ratio,rho_ratio] = isen_stag_ratios(M,gamma);
T = t_ratio.*T_0;
P = p_ratio.*p_0;

%%
res = 10;

subplot(2,1,1)
plot(x,r,'k-',x,r_e.*ones(size(x)),'k--',lip_coord(1),lip_coord(2),'kX');


axis equal
hold on

X_char = [];
Y_char = [];
T_char = [];
P_char = [];
for i = 1:length(x)
    x_cline = linspace(0,x(i),res);
    y_cline = (r(i)-lip_coord(2))./(x(i)-lip_coord(1)).*x_cline;
    T_cline = ones(size(x_cline)).*T(i);
    P_cline = ones(size(x_cline)).*P(i);
    X_char = [X_char; x_cline];
    Y_char = [Y_char; y_cline];
    T_char = [T_char; T_cline];
    P_char = [P_char; P_cline];
end

%contour(X_char',Y_char',T_char',100);
contour(X_char',Y_char',P_char',1000);

legend('Contour','Axisymmetric Centre','Nozzle Lip','Pressure Contour');

hold off
subplot(2,1,2)

hold on
plot(x,r,'k-',x,r_e.*ones(size(x)),'k--',lip_coord(1),lip_coord(2),'kX');

contour(X_char',Y_char',T_char',1000);

axis equal
legend('Contour','Axisymmetric Centre','Nozzle Lip','Temperature Contour');

hold off

csvwrite('aerospike_contour.csv',[x',M',T',P'])