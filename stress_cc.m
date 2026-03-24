% =========================================================================
% stress_cc.m
%
% Computes the radial stress component σ_χχ in bipolar coordinates (χ, ξ).
%
% This is the normal stress acting on surfaces of constant χ (i.e., on the
% circular hole boundaries and concentric circles). When evaluated at the
% hole boundary χ = χ₀, this should equal −P everywhere (boundary condition).
%
% The expression implements Eq. (A6) from Saeedi & Kothari (2025), which
% is a corrected version of the formula in Davanas (1992) and Ling (1948).
%
% Inputs:
%   chi  - bipolar radial coordinate (scalar or vector) [−]
%   xi   - bipolar angular coordinate (scalar or vector) [−]
%   A    - series coefficients A_n (1×N vector), from Eq. (A10)
%   B    - series coefficients B_n (1×N vector), from Eqs. (A11–A12)
%   K    - amplitude constant (scalar), from Eq. (A9)
%   N    - number of series terms (positive integer)
%   c    - bipolar focal distance (scalar) [same units as hole radius R]
%
% Output:
%   sig_cc - σ_χχ stress component at the given (chi, xi) coordinates
%            [same units as applied pressure P]
%
% Note: The final division by c converts from the non-dimensional series
% to the physical stress (see the c^{-1} factor in Eq. A6 of the paper).
% =========================================================================

function sig_cc = stress_cc(chi, xi, A, B, K, N, c)

    % Leading term: contribution from the constant K (non-series part of
    % the stress function). This captures the far-field behavior.
    sig_cc = -0.5 .* K .* (cosh(2.*chi) + cos(2.*xi) - 2.*cosh(chi).*cos(xi));

    % A_n series (n = 1 to N): each term contributes through hyperbolic
    % functions of (n+1)*chi, capturing the spatial decay away from the hole.
    for i = 1:N
        sig_cc = sig_cc + A(i) .* ( ...
            -(i)^2   .* cosh((i+1).*chi) .* (cosh(chi) - cos(xi)) .* cos(i.*xi) ...
            -(i+1)   .* sinh((i+1).*chi) .* sinh(chi)              .* cos(i.*xi) ...
            + i      .* sin(xi)          .* sin(i.*xi)              .* cosh((i+1).*chi) ...
            + cos(i.*xi) .* cosh(chi)    .* cosh((i+1).*chi) );
    end

    % B_1 term (constant): the n=1 B coefficient enters separately because
    % the B_n formula (Eq. A12) has a factor 1/(n*(n-1)) that is singular
    % at n=1; B_1 is given by the separate formula Eq. (A11).
    sig_cc = sig_cc + B(1);

    % B_n series (n = 2 to N): analogous to A_n but using (n-1)*chi,
    % capturing the interaction between the two holes.
    for i = 2:N
        sig_cc = sig_cc + B(i) .* ( ...
            -(i)^2   .* cosh((i-1).*chi) .* (cosh(chi) - cos(xi)) .* cos(i.*xi) ...
            -(i-1)   .* sinh((i-1).*chi) .* sinh(chi)              .* cos(i.*xi) ...
            + i      .* sin(xi)          .* sin(i.*xi)              .* cosh((i-1).*chi) ...
            + cos(i.*xi) .* cosh(chi)    .* cosh((i-1).*chi) );
    end

    % Normalize by the focal distance c to recover physical stress units
    sig_cc = sig_cc / c;

end
