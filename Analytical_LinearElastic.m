% Notation follows that of Davanas Paper with "a" replaced with "c" and
% "eta" replaced with "chi"
% equation of the hole on the right is c = R*sinh(chi1)
% The distance between the edges of the two holes is L = 2*(sqrt(c^2+R^2)-R)
% We want to keep the radius R constant, and vary the length L. This is
% only possible by changing c.

clear all
close all

% Material Parameters
R = 0.1;              % Radius of the hole
E = 52;             % Young's modulus
nu = 0.3;           % Poisson's ratio
P = 5;             % prescribed pressure
mu = E/2/(1+nu);     % Shear modulus

% Choose which case to run:
% case I: single value of chi at the internal surface
% case II: multiple values of chi (useful to plot the energy)

case_num = 2;

% Define a range of 'c' values (every 'c' corresponds to a hole)
if case_num==1
    half_eta = 0.3/2; % eta/2 = L/2 + R
    c_vec =  [sqrt(half_eta^2-R^2)];
elseif case_num==2
    c_vec = linspace(0.05,1 ,40);
else
    disp('invalid choice of case')
end

% create vectors for storing total strain energy and the edge-to-edge hole
% distance L
total_se = zeros(1, length(c_vec));
L = zeros(1,length(c_vec));

% preparation for figure
fig_stress = figure();
fig_stress.Position =  [0   0   560   1.5*560];
t = tiledlayout(3,1);
t.OuterPosition = [0 0 1 1]
ax1 = nexttile;
ax1.FontSize = 16;
ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.String = '$\theta$';

ax1.Title.String = 'Radial Stress';
ax1.LineWidth = 1.2;

ax2 = nexttile;
ax2.FontSize = 16;
ax2.XLabel.String = '\theta';
ax2.Title.String = 'Angular Stress';
ax2.LineWidth = 1.2;


ax3 = nexttile;
ax3.FontSize = 16;
ax3.XLabel.String = '\theta';
ax3.Title.String = 'Shear Stress';
ax3.LineWidth = 1.2;

h = gobjects(3,length(c_vec));

%%%%%%%%%%%%%%%%%% Calculation of energy for various values of distance(eta) by varying
%%%%%%%%%%%%%%%%%% the parameter c
for m=1:length(c_vec)

c = c_vec(m);

chi0 = asinh(c/R); % Corresponding chi value in the bipolar system

N = 100;

% defining the variables to store coefficients of the series A_n and B_n in
% Davanas, and in Ling 1948
A = zeros(1, N);
B = zeros(1,N);

S_cc = zeros(N,N);
S_xx = zeros(N,N);
S_cx = zeros(N,N);

%%% Vectors for the geometry to be used later

% For integration over entire space (all of xi and half of chi)
chi_vec = linspace(0, chi0,N);
xi_vec = linspace(-pi, pi, N);

% For checking stresses (xi from 0 to pi; that gives all the info because
% of symmetery in the problem)
xi_vec_plot = logspace(-10, pi, N);

% calculating K
F = 0; %initialize
for k = 2:1*N % just for higher accuracy; no other reason to go to 2N
    F = F + (exp(-k*chi0)*sinh(k*chi0) + k*sinh(chi0)*(k*sinh(chi0)+ cosh(chi0)))/(k*(k^2-1)*(sinh(2*k*chi0)+k*sinh(2*chi0)));
end

K = c*P*(0.5+ tanh(chi0)*(sinh(chi0))^2 - 4*F)^-1;

% Series coefficients calculation and verification
% calculating coeffiecients of stress function series in Davanas (using
% also help from Ling 1948 to correct the mistakes in both papers)
A_sum = 0;
for k = 1:N
    A(k) = 2*K*(exp(-k*chi0)*sinh(k*chi0) + k*exp(-chi0)*sinh(chi0))/(k*(k+1)*(sinh(2*k*chi0)+k*sinh(2*chi0)));
    A_sum  = A_sum + 2*K*(exp(-k*chi0)*sinh(k*chi0) + k*exp(-chi0)*sinh(chi0))/(k*(k+1)*(sinh(2*k*chi0)+k*sinh(2*chi0)));
end

B_sum = 0.5*(K*tanh(chi0)*cosh(2*chi0) - 2*c*P);
B(1) =  0.5*(K*tanh(chi0)*cosh(2*chi0) - 2*c*P);
for k = 2:N
    B(k) = - 2*K*(exp(-k*chi0)*sinh(k*chi0) + k*exp(chi0)*sinh(chi0))/(k*(k-1)*(sinh(2*k*chi0)+k*sinh(2*chi0)));
    B_sum = B_sum - 2*K*(exp(-k*chi0)*sinh(k*chi0) + k*exp(chi0)*sinh(chi0))/(k*(k-1)*(sinh(2*k*chi0)+k*sinh(2*chi0)));
end

% Verification#1: Verifying that A_sum+B_sum should be zero.
if(abs(A_sum + B_sum) <1E-6)
    disp('Passed: Sum of coefficients passed the check')
else
    disp('Failed: Sum of coefficients failed the check')
end

% Verification#2: Verifying that stress_cc at hole = -P and stress_cx at hole = 0
s_radial_boundary = stress_cc(chi0, xi_vec_plot, A, B, K, N, c);
s_shear_boundary = stress_cx(chi0, xi_vec_plot, A, B, K, N, c);

if ((sum(s_radial_boundary.*s_radial_boundary) - P^2*N <1E-6) && (sum(s_shear_boundary.*s_shear_boundary) <1E-6))
    disp('Passed: Stresses are zero at the internal boundary')
else
    disp('Failed: Stresses are not zero at the internal boundary')
    s_radial_boundary;
    s_shear_boundary;
end

%%% Calculating strain energy

%calculation 1 / brute force and works
se_mat = zeros(N,N);
for i = 1:N
    for j = 1:N
        se_mat(i,j) = (stress_cc(chi_vec(i),xi_vec(j),A,B,K,N,c).*(stress_cc(chi_vec(i),xi_vec(j),A,B,K,N,c)/E*(1-nu^2) - (1+nu)*nu*(stress_xx(chi_vec(i),xi_vec(j),A,B,K,N,c))/E)+...
            stress_xx(chi_vec(i),xi_vec(j),A,B,K,N,c).*(stress_xx(chi_vec(i),xi_vec(j),A,B,K,N,c)/E*(1-nu^2) - (1+nu)*nu*(stress_cc(chi_vec(i),xi_vec(j),A,B,K,N,c))/E)+...
            stress_cx(chi_vec(i),xi_vec(j),A,B,K,N,c).*stress_cx(chi_vec(i),xi_vec(j),A,B,K,N,c)/E*2*(1+nu)).*(c./(cosh(chi_vec(i))-cos(xi_vec(j)))).^2;
        
        % for plotting stress later
        S_cc(i,j) = stress_cc(chi_vec(i),xi_vec(j),A,B,K,N,c);
        S_xx(i,j) = stress_xx(chi_vec(i),xi_vec(j),A,B,K,N,c);
        S_cx(i,j) = stress_cx(chi_vec(i),xi_vec(j),A,B,K,N,c);
    end
end

% %calculation 2 /vectorization
se_mat2 = zeros(N,N);
for j = 1:N
        se_mat2(:,j) = (stress_cc(chi_vec,xi_vec(j),A,B,K,N,c).*(stress_cc(chi_vec,xi_vec(j),A,B,K,N,c)/E*(1-nu^2) - (1+nu)*nu*(stress_xx(chi_vec,xi_vec(j),A,B,K,N,c))/E)+...
            stress_xx(chi_vec,xi_vec(j),A,B,K,N,c).*(stress_xx(chi_vec,xi_vec(j),A,B,K,N,c)/E*(1-nu^2) - (1+nu)*nu*(stress_cc(chi_vec,xi_vec(j),A,B,K,N,c))/E)+...
            stress_cx(chi_vec,xi_vec(j),A,B,K,N,c).*stress_cx(chi_vec,xi_vec(j),A,B,K,N,c)/E*2*(1+nu)).*c./(cosh(chi_vec)-cos(xi_vec(j)));
end

total_se(m) =  trapz(chi_vec,trapz(xi_vec, se_mat,2));
total_se2(m) =  trapz(chi_vec,trapz(xi_vec, se_mat2,2));
L(m) = 2*(sqrt(c^2+R^2)-R);

%%% Calculating the stress at another circle given by chi = constant
chi_test = 0.9*chi0; % chi0 is the hole. All circles outside the hole are represented by chi<chi0.

% Convert from xi to theta (they are not the same even though they both go
% from [0,pi]. theta changes non-linearly as a function of xi.
theta_vec = acos((-1+cosh(chi_test)*cos(xi_vec_plot))./(cosh(chi_test)-cos(xi_vec_plot)));

% Calculating the stresses at the test case
s_radial_test = stress_cc(chi_test, theta_vec, A, B, K, N, c);
s_angular_test = stress_xx(chi_test, theta_vec , A, B, K, N, c);
s_shear_test = stress_cx(chi_test, theta_vec, A, B, K, N, c);

%% Plotting stress every iteration
axes(ax1);
hold(ax1,"on")
h(1,m) = plot(ax1, theta_vec, s_radial_test, 'LineWidth', 1.5);
hold(ax1,"off")

axes(ax2);
hold(ax2,"on")
h(2,m) = plot(ax2, theta_vec, s_angular_test, 'LineWidth', 1.5);
hold(ax2,"off")

axes(ax3);
hold(ax3,"on")
h(3,m) = plot(ax3, theta_vec, s_shear_test, 'LineWidth', 1.5);
hold(ax3,"off")
end

lgd = legend(h(1,:),[num2str(c_vec')]);
title(lgd, 'c=');
lgd.Layout.Tile = 'East'; % to place the legend in tiled system, this is the property to use.

energy_inf = pi*R^2*P^2/mu;
%%% Plotting

% Total energy as a function of center-to-center hole distance
% SE vs eta/R = (L+2R)/R
figure;
plot(L/R +2, total_se/energy_inf, 'om', 'LineWidth',2, 'MarkerSize', 5)
hold on
plot(L/R +2, total_se/energy_inf, '*k', 'LineWidth',2, 'MarkerSize', 5)
grid on
xlabel('Distance to radius ratio ($\eta$/R)','Interpreter','latex')
ylabel('Nondim. SE ($SE_{total}/{|SE_\infty|}$)','Interpreter','latex')
legend('First Approach - Brute Force','Second Approach - Vectorization')

% Calculate analytical PE
dR_anlyt = L/R +2;
PE = -total_se;
PE_inf = -energy_inf;
PEn_analyt = PE/abs(PE_inf);
figure
plot(dR_anlyt, PEn_analyt, 'ok', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Distance to radius ratio ($\eta$/R)','Interpreter','latex');
ylabel('Nondim. PE ($PE_{total}/{|PE_\infty|}$)','Interpreter','latex');
grid on

% Calculate first derivative using numerical differentiation
figure
plot(dR_anlyt, -gradient(PEn_analyt, dR_anlyt)*abs(PE_inf)/mu/R, 'ok', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Distance to radius ratio ($\eta$/R)','Interpreter','latex');
ylabel('Nondim. Driving Force ($\mathcal{F}/\mu$)','Interpreter','latex');
grid on;



