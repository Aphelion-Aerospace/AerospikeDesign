function [X,Y,theta,nu,M,mu,K_neg,K_pos] = anderson_nozzle_contour(A_t,epsilon, gamma,points)
r_t = sqrt(A_t./pi);
n = 7;

X = zeros(1/2.*n.*(n+3),1);
Y = r_t.*ones(size(X));
theta = zeros(size(X));
nu = zeros(size(X));
M = zeros(size(X));
mu = zeros(size(X));
K_neg = zeros(size(X));
K_pos = zeros(size(X));



A_e = A_t.*epsilon;

[M_e,f_dummy] = fsolve(@(M) mach_expansion_ratio(M,gamma,epsilon),10);

M_e =2.4
%% Initial expansion waves
theta_init = linspace(0.375,18.375,n);

theta(1:n) = theta_init.*pi./180;

M_init = zeros(size(theta_init));

for i = 1:length(theta_init)
    [M_init(i), func_dummy] = fsolve(@(M) prandtl_meyer_zero(M,gamma,theta_init(i)),10);
end

M_init = real(M_init);

M(1:length(M_init)) = M_init;
nu(1:length(theta_init)) = prandtl_meyer(M_init,gamma);

mu(1:length(theta_init)) = mach_angle(M_init);

K_neg(1:length(theta_init)) = theta(1:length(theta_init))+nu(1:length(theta_init));
K_pos(1:length(theta_init)) = theta(1:length(theta_init))-nu(1:length(theta_init));


theta(n+1) = theta(n);
nu(n+1) = nu(n);
M(n+1) = M(n);
mu(n+1) = mu(n);
K_neg(n+1) = K_neg(n);
K_pos(n+1) = K_pos(n);

[X(1),Y(1)] = line_int(0,0,0,theta(1),0,r_t);

for i = 2:n
    [X(i),Y(i)] = line_int(theta(i),0,r_t,1./2.*(theta(i-1)+theta(i)) + 1./2.*(mu(i-1)+mu(i)),X(i-1),Y(i-1))
end
[X(n+1),Y(n+1)] = line_int(1./2.*(18.375 + theta(n+1)),X(1),Y(1),1./2.*(theta(n+1)+theta(n))+1./2.*(mu(n+1)+mu(n)),X(n),Y(n))
%%
%centre line points
cPts = ones(n,1);
count = 1;
i = 1
while (i < 1./2.*n.*(n+3))
    cPts(count) = i;
    count = count + 1;
    i = i + n - count + 3;
end

%nozzle pts
nPts = ones(n,1);
count = 1;
i = 1;
while (count < n)
    nPts(count) = n + i - count + 1
    count = count + 1;
    i = cPts(count);
end
nPts(end) = 1./2.*n.*(n+3);

cPts = cPts(2:end);
nPts = nPts(2:end);

for i = (n+1):(1./2.*n.*(n+3))
    ajm = n - ceil(i/n)-2;
    if ismember(i,nPts)
        theta(i) = theta(i-1);
        nu(i) = nu(i-1);
        M(i) = M(i-1);
        mu(i) = mu(i-1);
        K_neg(i) = K_neg(i-1);
        K_pos(i) = K_pos(i-1);
        [X(i),Y(i)] = line_int(1./2.*(theta(i)+theta(i-ajm)),X(i-ajm),Y(i-ajm),1./2.*(theta(i)+theta(i-1))+1./2.*(mu(i)+mu(i-1)),X(i-1),Y(i-1));
    elseif ismember(i,cPts)
        theta(i) = 0;
        nu(i) = nu(i-1)+theta(i-1);
        [M(i), f_dummy] = fsolve(@(M) prandtl_meyer_zero(M,gamma,nu(i)),10);
        mu(i) = mach_angle(M(i));
        K_neg(i) = theta(i) + nu(i);
        K_pos(i) = theta(i) - nu(i);
        [X(i),Y(i)] = line_int(0,0,0,1./2.*(theta(i) + theta(i-ajm)) - 1./2.*(mu(i)+ mu(i-ajm)),X(i-ajm),Y(i-ajm));
    else
        theta(i) = 1./2.*(K_neg(i-ajm) + K_pos(i-1));
        nu(i) = 1./2.*(K_neg(i-ajm) - K_pos(i-1));
        [M(i), f_dummy] = fsolve(@(M) prandtl_meyer_zero(M,gamma,nu(i)),10);
        mu(i) = mach_angle(M(i));
        K_neg(i) = theta(i) + nu(i);
        K_pos(i) = theta(i) - nu(i);
        [X(i),Y(i)] = line_int(1./2.*(theta(i-ajm) + theta(i)) - 1./2.*(mu(i) + mu(i-ajm)),X(i-ajm),Y(i-ajm),1./2.*(theta(i) +theta(i-1))+1./2.*(mu(i)+mu(i-1)),X(i-1),Y(i-1));
    end
end
