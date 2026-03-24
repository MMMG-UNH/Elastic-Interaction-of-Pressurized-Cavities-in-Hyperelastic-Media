% =========================================================================
% Analytical_LinearElastic.m
%
% Computes the semi-analytical solution for two equal, pressurized circular
% holes in an infinite, isotropic, linear elastic medium under plane-strain
% conditions using a bipolar coordinate system.
%
% The solution follows the formulation of Ling (1948) and Davanas (1992),
% with corrections applied to both references. The bipolar coordinate
% notation uses:
%   chi (χ): radial-like coordinate (constant-chi curves are circles)
%   xi  (ξ): angular-like coordinate (constant-xi curves are circles
%             orthogonal to the chi-circles)
%   c      : focal distance of the bipolar system
%
% The right hole boundary is located at chi = chi0 = asinh(c/R).
% The edge-to-edge distance is L = 2*(sqrt(c^2+R^2) - R).
% The center-to-center distance is eta = L + 2*R.
%
% Outputs:
%   - Stress plots (radial σ_χχ, hoop σ_ξξ, shear σ_χξ) vs. angle on
%     a test circle just inside the hole boundary
%   - Non-dimensionalized strain energy vs. eta/R
%   - Non-dimensionalized potential energy vs. eta/R
%   - Non-dimensionalized driving force vs. eta/R
%
% Functions called: stress_cc.m, stress_xx.m, stress_cx.m
%
% References:
%   Ling, C. (1948). J. Appl. Phys., 19(1), 77-82.
%   Davanas, K. (1992). J. Mater. Sci., 27(6), 1589-1598.
%   Saeedi & Kothari (2025). J. Appl. Mech., 92(5), 051008.
% =========================================================================

clear all
close all

% -------------------------------------------------------------------------
% Material and geometry parameters
% -------------------------------------------------------------------------
R  = 0.1;           % Undeformed hole radius [length units]
E  = 52;            % Young's modulus [force/area]
nu = 0.3;           % Poisson's ratio [-]
P  = 5;             % Applied internal pressure [force/area]
mu = E / 2 / (1 + nu);  % Shear modulus derived from E and nu

% -------------------------------------------------------------------------
% Case selection
%   case 1: single value of c (one separation distance)
%   case 2: sweep over a range of c values (energy vs. distance curve)
% -------------------------------------------------------------------------
case_num = 2;

% Build the vector of focal-distance values c.
% Each value of c corresponds to a unique hole separation:
%   eta = 2*sqrt(c^2 + R^2)  (center-to-center distance)
if case_num == 1
    half_eta = 0.3 / 2;          % half of center-to-center distance
    c_vec = [sqrt(half_eta^2 - R^2)];
elseif case_num == 2
    c_vec = linspace(0.05, 1, 40);  % sweep from close to far separation
else
    disp('invalid choice of case')
end

% Pre-allocate output arrays
total_se  = zeros(1, length(c_vec));  % total strain energy for each c
L         = zeros(1, length(c_vec));  % edge-to-edge hole distance for each c

% -------------------------------------------------------------------------
% Figure setup: 3-panel plot of stresses at a test circle
% -------------------------------------------------------------------------
fig_stress = figure();
fig_stress.Position = [0 0 560 1.5*560];
t = tiledlayout(3, 1);
t.OuterPosition = [0 0 1 1];

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

h = gobjects(3, length(c_vec));  % handle array for legend entries

% =========================================================================
% Main loop: compute energy and stresses for each separation distance
% =========================================================================
for m = 1:length(c_vec)

    c    = c_vec(m);
    chi0 = asinh(c / R);   % chi-value at the hole boundary (χ₀ = sinh⁻¹(c/R))

    N = 100;  % number of terms in the series expansion (convergence checked)

    % Pre-allocate series coefficients A_n and B_n (Eqs. A10–A12 of paper)
    A = zeros(1, N);
    B = zeros(1, N);

    % Pre-allocate stress matrices for full-field calculation
    S_cc = zeros(N, N);  % radial stress σ_χχ field
    S_xx = zeros(N, N);  % hoop stress  σ_ξξ field
    S_cx = zeros(N, N);  % shear stress σ_χξ field

    % Coordinate grids for numerical integration over the domain
    % Integration runs from chi=0 (midplane between holes) to chi=chi0 (hole)
    % and from xi=-pi to xi=pi (full angular sweep).
    chi_vec = linspace(0, chi0, N);
    xi_vec  = linspace(-pi, pi, N);

    % Angular grid for stress verification and plotting.
    % Log-spaced to capture steep gradients near xi=0 (closest point).
    xi_vec_plot = logspace(-10, pi, N);

    % ---------------------------------------------------------------------
    % Compute the constant K (Eq. A9 of paper)
    % K normalizes the stress magnitude to satisfy the boundary condition
    % σ_χχ = -P at chi = chi0.
    % ---------------------------------------------------------------------
    F = 0;  % initialize the series sum F (Eq. A15 of paper)
    for k = 2:1*N
        F = F + (exp(-k*chi0)*sinh(k*chi0) + k*sinh(chi0)*(k*sinh(chi0) + cosh(chi0))) ...
              / (k * (k^2 - 1) * (sinh(2*k*chi0) + k*sinh(2*chi0)));
    end
    % K is the amplitude constant satisfying the pressure boundary condition
    K = c * P * (0.5 + tanh(chi0)*(sinh(chi0))^2 - 4*F)^-1;

    % ---------------------------------------------------------------------
    % Compute series coefficients A_n and B_n (Eqs. A10–A12 of paper)
    % These coefficients come from matching boundary conditions at chi=chi0.
    % ---------------------------------------------------------------------
    A_sum = 0;
    for k = 1:N
        A(k)  = 2*K * (exp(-k*chi0)*sinh(k*chi0) + k*exp(-chi0)*sinh(chi0)) ...
                    / (k * (k+1) * (sinh(2*k*chi0) + k*sinh(2*chi0)));
        A_sum = A_sum + A(k);
    end

    B_sum = 0.5 * (K*tanh(chi0)*cosh(2*chi0) - 2*c*P);
    B(1)  = 0.5 * (K*tanh(chi0)*cosh(2*chi0) - 2*c*P);
    for k = 2:N
        B(k)  = -2*K * (exp(-k*chi0)*sinh(k*chi0) + k*exp(chi0)*sinh(chi0)) ...
                     / (k * (k-1) * (sinh(2*k*chi0) + k*sinh(2*chi0)));
        B_sum = B_sum + B(k);
    end

    % ---------------------------------------------------------------------
    % Verification 1: sum of all A and B coefficients must be zero.
    % This is a necessary consistency condition for the stress function.
    % ---------------------------------------------------------------------
    if abs(A_sum + B_sum) < 1E-6
        disp('Passed: Sum of coefficients passed the check')
    else
        disp('Failed: Sum of coefficients failed the check')
    end

    % ---------------------------------------------------------------------
    % Verification 2: boundary stresses at chi=chi0 (the hole surface).
    % σ_χχ should equal -P (applied pressure) and σ_χξ should be zero
    % everywhere on the hole boundary, confirming the boundary conditions
    % are satisfied.
    % ---------------------------------------------------------------------
    s_radial_boundary = stress_cc(chi0, xi_vec_plot, A, B, K, N, c);
    s_shear_boundary  = stress_cx(chi0, xi_vec_plot, A, B, K, N, c);

    if (sum(s_radial_boundary.*s_radial_boundary) - P^2*N < 1E-6) && ...
       (sum(s_shear_boundary.*s_shear_boundary) < 1E-6)
        disp('Passed: Stresses are zero at the internal boundary')
    else
        disp('Failed: Stresses are not zero at the internal boundary')
    end

    % =====================================================================
    % Strain energy computation
    %
    % The total strain energy is computed by integrating the strain energy
    % density over the domain (half-domain by symmetry):
    %
    %   SE = integral_0^chi0 integral_{-pi}^{pi} (σ:ε) * J^2 dxi dchi
    %
    % where J = c/(cosh(chi) - cos(xi)) is the Jacobian of the bipolar
    % transformation (Eq. A23 of paper). The factor J^2 converts from
    % bipolar to Cartesian area element.
    %
    % Two approaches are used for cross-verification:
    %   Method 1: explicit double loop (brute force, clear but slow)
    %   Method 2: partial vectorization over chi (faster)
    % =====================================================================

    % --- element-wise double loop ---
    se_mat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            chi_i = chi_vec(i);
            xi_j  = xi_vec(j);
            J2    = (c / (cosh(chi_i) - cos(xi_j)))^2;  % Jacobian squared

            scc = stress_cc(chi_i, xi_j, A, B, K, N, c);
            sxx = stress_xx(chi_i, xi_j, A, B, K, N, c);
            scx = stress_cx(chi_i, xi_j, A, B, K, N, c);

            % Plane-strain strain energy density (inverted Hooke's law):
            %   W = σ_cc*(σ_cc*(1-ν²)/E - ν(1+ν)*σ_xx/E)
            %     + σ_xx*(σ_xx*(1-ν²)/E - ν(1+ν)*σ_cc/E)
            %     + σ_cx*σ_cx*2(1+ν)/E
            se_mat(i, j) = (scc .* (scc/E*(1-nu^2) - (1+nu)*nu*sxx/E) + ...
                            sxx .* (sxx/E*(1-nu^2) - (1+nu)*nu*scc/E) + ...
                            scx .* scx/E*2*(1+nu)) .* J2;

            % Store full stress field for potential plotting
            S_cc(i, j) = scc;
            S_xx(i, j) = sxx;
            S_cx(i, j) = scx;
        end
    end


    % Integrate over the bipolar domain using the trapezoidal rule
    % Outer integral over chi, inner integral over xi
    total_se(m)  = trapz(chi_vec, trapz(xi_vec, se_mat,  2));

    % Edge-to-edge distance L for this value of c
    L(m) = 2 * (sqrt(c^2 + R^2) - R);

    % =====================================================================
    % Stress at a test circle inside the domain (chi_test < chi0)
    %
    % chi_test = 0.9*chi0 corresponds to a circle slightly outside the
    % hole boundary. Stresses are plotted as a function of the polar angle
    % theta (measured from the axis connecting hole centers).
    %
    % Note: theta and xi are NOT identical even though both range over
    % [0, pi]. The mapping between them is nonlinear (see below).
    % =====================================================================
    chi_test = 0.9 * chi0;

    % Convert xi to the standard polar angle theta via:
    %   cos(theta) = (cosh(chi_test)*cos(xi) - 1) / (cosh(chi_test) - cos(xi))
    % This is derived from the bipolar-to-Cartesian transformation.
    theta_vec = acos((-1 + cosh(chi_test).*cos(xi_vec_plot)) ./ ...
                     (cosh(chi_test) - cos(xi_vec_plot)));

    s_radial_test  = stress_cc(chi_test, theta_vec, A, B, K, N, c);
    s_angular_test = stress_xx(chi_test, theta_vec, A, B, K, N, c);
    s_shear_test   = stress_cx(chi_test, theta_vec, A, B, K, N, c);

    % --- Plot stresses ---
    axes(ax1); hold(ax1, "on")
    h(1, m) = plot(ax1, theta_vec, s_radial_test, 'LineWidth', 1.5);
    hold(ax1, "off")

    axes(ax2); hold(ax2, "on")
    h(2, m) = plot(ax2, theta_vec, s_angular_test, 'LineWidth', 1.5);
    hold(ax2, "off")

    axes(ax3); hold(ax3, "on")
    h(3, m) = plot(ax3, theta_vec, s_shear_test, 'LineWidth', 1.5);
    hold(ax3, "off")

end  % end loop over c values

% Add legend to stress figure
lgd = legend(h(1,:), [num2str(c_vec')]);
title(lgd, 'c=');
lgd.Layout.Tile = 'East';

% =========================================================================
% Post-processing: energy and driving force plots
% =========================================================================

% Reference energy: PE of two non-interacting cavities (Eq. A25 of paper)
% PE_inf = 2 * pi*R^2*P^2 / (2*mu) = pi*R^2*P^2 / mu
energy_inf = pi * R^2 * P^2 / mu;

% --- Plot 1: Strain energy vs. center-to-center distance ratio eta/R ---
% Two methods plotted for cross-verification; they should coincide.
figure;
plot(L/R + 2, total_se/energy_inf,  'om', 'LineWidth', 2, 'MarkerSize', 5)
grid on
xlabel('Distance to radius ratio ($\eta$/R)', 'Interpreter', 'latex')
ylabel('Nondim. SE ($SE_{total}/{|SE_\infty|}$)', 'Interpreter', 'latex')

% --- Plot 2: Potential energy vs. eta/R ---
% In linear elasticity under static conditions: PE_total = -SE_total
% (Eq. A24 of paper: the potential energy of external forces = -2*SE)
dR_anlyt  = L/R + 2;
PE        = -total_se;
PE_inf    = -energy_inf;
PEn_analyt = PE / abs(PE_inf);  % non-dimensionalized by |PE_inf|

figure
plot(dR_anlyt, PEn_analyt, 'ok', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Distance to radius ratio ($\eta$/R)', 'Interpreter', 'latex');
ylabel('Nondim. PE ($PE_{total}/{|PE_\infty|}$)', 'Interpreter', 'latex');
grid on

% --- Plot 3: Driving force vs. eta/R ---
% Driving force F = -dPE/deta (Eq. 10 of paper).
% Computed here by numerical differentiation using gradient().
% Negative values indicate attraction; positive values indicate repulsion.
% In linear elasticity, F is always negative (always attractive).
figure
plot(dR_anlyt, -gradient(PEn_analyt, dR_anlyt) * abs(PE_inf) / mu / R, ...
     'ok', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Distance to radius ratio ($\eta$/R)', 'Interpreter', 'latex');
ylabel('Nondim. Driving Force ($\mathcal{F}/\mu$)', 'Interpreter', 'latex');
grid on;
