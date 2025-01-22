% Notation follows that of Davanas Paper
% equation of the hole on the right is a = R*sinh(eta1)
% The distance between the edges of the two holes is L = 2*(sqrt(a^2+R^2)-R)
% We want to keep the radius R constant, and vary the length L. This is
% only possible by changing a.

clear all
%close all

% Material Parameters
R = 0.1;              % Radius of the hole
E = 52;             % Young's modulus
nu = 0.3;           % Poisson's ratio
P0 = 5;             % prescribed pressure
G = E/2/(1+nu);     % Shear modulus
eta_test = 1;       % choose the circle where the stress is evaluated


% Choose which case to run:
% case I: single value of eta at the internal surface
% case II: multiple values of eta (useful to plot the energy)

case_num = 2;


% Define a range of 'a' values (every 'a' corresponds to a hole)
if case_num==1
    a_vec =  [sqrt(2.5^2-1)];
elseif case_num ==2
    % a_vec  = [sqrt(2.5^2-1) sqrt(3.5^2-1)]
    % a_vec = linspace(R*sinh(eta_test),8 ,10);
    a_vec = linspace(0.05,1 ,40);
else
    disp('invalid choice of case')
end

% create vectors for storing total strain energy and the edge-to-edge hole
% distance L
total_se = zeros(1, length(a_vec));
L = zeros(1,length(a_vec));


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

h = gobjects(3,length(a_vec));

%%%%%%%%%%%%%%%%%% Calculation of energy for various values of d by varying
%%%%%%%%%%%%%%%%%% the parameter a
for m=1:length(a_vec)

a = a_vec(m);

eta0 = asinh(a/R); % Corresponding eta value in the bipolar system

N = 100;

% defining the variables to store coefficients of the series A_n and B_n in
% Davanas, and in Ling 1948
A = zeros(1, N);
B = zeros(1,N);

x = zeros(N,N);
S_ee = zeros(N,N);
S_xx = zeros(N,N);
S_ex = zeros(N,N);

%%% Vectors for the geometry to be used later

% For integration over entire space (all of xi and half of eta)
eta_vec = linspace(0, eta0,N);
xi_vec = linspace(-pi, pi, N);

% For checking stresses (xi from 0 to pi; that gives all the info because
% of symmetery in the problem)
xi_vec_plot = logspace(-10, pi, N);

% calculating K
F = 0; %initialize
for k = 2:1*N % just for higher accuracy; no other reason to go to 2N
    F = F+ (exp(-k*eta0)*sinh(k*eta0) + k*sinh(eta0)*(k*sinh(eta0)+ cosh(eta0)))/(k*(k^2-1)*(sinh(2*k*eta0)+k*sinh(2*eta0)));
end

K = a*P0*(0.5+ tanh(eta0)*(sinh(eta0))^2 - 4*F)^-1;


% Series coefficients calculation and verification
% calculating coeffiecients of stress function series in Davanas (using
% also help from Ling 1948 to correct the mistakes in both papers)
A_sum = 0;
for k = 1:N
    A(k) = 2*K*(exp(-k*eta0)*sinh(k*eta0) + k*exp(-eta0)*sinh(eta0))/(k*(k+1)*(sinh(2*k*eta0)+k*sinh(2*eta0)));
    A_sum  = A_sum + 2*K*(exp(-k*eta0)*sinh(k*eta0) + k*exp(-eta0)*sinh(eta0))/(k*(k+1)*(sinh(2*k*eta0)+k*sinh(2*eta0)));
end

B_sum = 0.5*(K*tanh(eta0)*cosh(2*eta0) - 2*a*P0);
B(1) =  0.5*(K*tanh(eta0)*cosh(2*eta0) - 2*a*P0);
for k = 2:N
    B(k) = - 2*K*(exp(-k*eta0)*sinh(k*eta0) + k*exp(eta0)*sinh(eta0))/(k*(k-1)*(sinh(2*k*eta0)+k*sinh(2*eta0)));
    B_sum = B_sum - 2*K*(exp(-k*eta0)*sinh(k*eta0) + k*exp(eta0)*sinh(eta0))/(k*(k-1)*(sinh(2*k*eta0)+k*sinh(2*eta0)));
end

% Verification#1: Verifying that A_sum+B_sum should be zero.
if(abs(A_sum + B_sum) <1E-6)
    disp('Passed: Sum of coefficients passed the check')
else
    disp('Failed: Sum of coefficients failed the check')
end

% Verification#2: Verifying that stress_ee at hole = -P0 and stress_ex at hole = 0
s_radial_boundary = stress_ee(eta0, xi_vec_plot, A, B, K, N,a);
s_shear_boundary=stress_ex(eta0, xi_vec_plot, A, B, K, N,a);

if ((sum(s_radial_boundary.*s_radial_boundary)- P0^2*N <1E-6) && (sum(s_shear_boundary.*s_shear_boundary) <1E-6))
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
        se_mat(i,j) = (stress_ee(eta_vec(i),xi_vec(j),A,B,K,N,a).*(stress_ee(eta_vec(i),xi_vec(j),A,B,K,N,a)/E*(1-nu^2) - (1+nu)*nu*(stress_xx(eta_vec(i),xi_vec(j),A,B,K,N,a))/E)+...
            stress_xx(eta_vec(i),xi_vec(j),A,B,K,N,a).*(stress_xx(eta_vec(i),xi_vec(j),A,B,K,N,a)/E*(1-nu^2) - (1+nu)*nu*(stress_ee(eta_vec(i),xi_vec(j),A,B,K,N,a))/E)+...
            stress_ex(eta_vec(i),xi_vec(j),A,B,K,N,a).*stress_ex(eta_vec(i),xi_vec(j),A,B,K,N,a)/E*2*(1+nu)).*(a./(cosh(eta_vec(i))-cos(xi_vec(j)))).^2;
        
        % for plotting stress later
        x(i,j) = a./(cosh(eta_vec(i))-cos(xi_vec(j)))*sinh(eta_vec(i));
        S_ee(i,j) = stress_ee(eta_vec(i),xi_vec(j),A,B,K,N,a);
        S_xx(i,j) = stress_xx(eta_vec(i),xi_vec(j),A,B,K,N,a);
        S_ex(i,j) = stress_ex(eta_vec(i),xi_vec(j),A,B,K,N,a);
    end
end

% %calculation 2 /vectorization
se_mat2 = zeros(N,N);
for j = 1:N
        se_mat2(:,j) = (stress_ee(eta_vec,xi_vec(j),A,B,K,N,a).*(stress_ee(eta_vec,xi_vec(j),A,B,K,N,a)/E*(1-nu^2) - (1+nu)*nu*(stress_xx(eta_vec,xi_vec(j),A,B,K,N,a))/E)+...
            stress_xx(eta_vec,xi_vec(j),A,B,K,N,a).*(stress_xx(eta_vec,xi_vec(j),A,B,K,N,a)/E*(1-nu^2) - (1+nu)*nu*(stress_ee(eta_vec,xi_vec(j),A,B,K,N,a))/E)+...
            stress_ex(eta_vec,xi_vec(j),A,B,K,N,a).*stress_ex(eta_vec,xi_vec(j),A,B,K,N,a)/E*2*(1+nu)).*a./(cosh(eta_vec)-cos(xi_vec(j)));
end


total_se(m) =  trapz(eta_vec,trapz(xi_vec, se_mat,2));
total_se2(m) =  trapz(eta_vec,trapz(xi_vec, se_mat2,2));
L(m) = 2*(sqrt(a^2+R^2)-R);


%%% Calculating the stress at another circle given by eta = constant
eta_test = 0.9*eta0; % eta0 is the hole. All circles outside the hole are represented by eta<eta0.


% Convert from xi to theta (they are not the same even though they both go
% from [0,pi]. Theta changes non-linearly as a function of xi.
% theta_vec = asin(sin(xi_vec_plot).*sinh(eta_test)./(cosh(eta_test)-cos(xi_vec_plot)));

theta_vec = acos((-1+cosh(eta_test)*cos(xi_vec_plot))./(cosh(eta_test)-cos(xi_vec_plot)));


% Calculating the stresses at the test case
s_radial_test = stress_ee(eta_test, theta_vec, A, B, K, N,a);
s_angular_test =stress_xx(eta_test, theta_vec , A, B, K, N,a);
s_shear_test =stress_ex(eta_test, theta_vec, A, B, K, N,a);


%% Plotting stress every iteration
axes(ax1);
hold(ax1,"on")
h(1,m) = plot(ax1, theta_vec, s_radial_test, 'LineWidth', 1.5);
hold(ax1,"off")

axes(ax2);
hold(ax2,"on")
h(2,m) = plot(ax2, theta_vec,s_angular_test, 'LineWidth', 1.5);
hold(ax2,"off")

axes(ax3);
hold(ax3,"on")
h(3,m)=plot(ax3,theta_vec,s_shear_test, 'LineWidth', 1.5);
hold(ax3,"off")
end

lgd = legend(h(1,:),[num2str(a_vec')]);
title(lgd, 'a=');
lgd.Layout.Tile = 'East'; % to place the legend in tiled system, this is the property to use.


energy_inf = pi*R^2*P0^2/G;
%%% Plotting

% Total energy as a function of center-to-center hole distance
% SE vs d
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
plot(dR_anlyt, -gradient(PEn_analyt, dR_anlyt)*abs(PE_inf)/G/R, 'ok', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Distance to radius ratio ($\eta$/R)','Interpreter','latex');
ylabel('Nondim. Driving Force ($\mathcal{F}/\mu$)','Interpreter','latex');
grid on;



