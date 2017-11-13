% Design aerospike nozzle for 20k

% Design params.
r_e = 0.0762;
epsilon = 8.1273;
A_t = r_e.^2.*pi./epsilon; % max expansion (r_b = 0) (r_e.^2 >= A_t.*epsilon./pi)

% Taken from CEA
gamma = mean([1.2534,1.2852]); % mean of throat and exit gamma
T_0 = 2833.63; % K
P_0 = 34.474; % bar
cp = mean([1.8292,1.6196]);
R = cp.*(1-1./gamma);

% Getting discrete contour points
[x,r,M,A] = angelino_nozzle_contour(A_t,epsilon,r_e,gamma,1000); %100 approximation points
lip_coord = [0,0];

% Calculating flow properties
[t_ratio,p_ratio,rho_ratio] = isen_stag_ratios(M,gamma);
T = T_0.*t_ratio; P = P_0.*p_ratio;
a = sqrt(gamma.*R.*T);
V = M.*a;
% Transforming to cylindrical coord
r = r_e-r;
lip_coord(2) = r_e;

lip_coord(1) = lip_coord(1) - min(x);
x = x - min(x);

plot(x,r,lip_coord(1),lip_coord(2),'o')
axis equal

% Getting line integral form

[x,s] = convert_to_path_length(x,r);

%% Outputting results to CSV

csvwrite('aerospike_contour.csv',[x',r',s',P',T',M',A',a',V'])
