% =========================================================================
% stress_xx.m
%
% Computes the hoop stress component σ_ξξ in bipolar coordinates (χ, ξ).
%
% This is the normal stress acting on surfaces of constant ξ (orthogonal
% circles in the bipolar system). Together with σ_χχ and σ_χξ, it fully
% characterizes the in-plane stress state.
%
% The expression implements Eq. (A7) from Saeedi & Kothari (2025), which
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
%   sig_xx - σ_ξξ stress component at the given (chi, xi) coordinates
%            [same units as applied pressure P]
%
% Note: The key structural difference from stress_cc.m is:
%   - The leading K term has opposite sign (+½ vs −½)
%   - The A_n/B_n terms use (n+1)² and cos(ξ)*cosh instead of cosh(χ)*cosh
%   These differences reflect the distinct forms of σ_χχ vs σ_ξξ in the
%   Airy stress function framework (compare Eqs. A6 and A7).
% =========================================================================

function sig_xx = stress_xx(chi, xi, A, B, K, N, c)

    % Leading term: note the positive sign (+½K), opposite to stress_cc.
    % This reflects that σ_χχ + σ_ξξ = K*(something), which is a consequence
    % of the biharmonic Airy stress function formulation.
    sig_xx = 0.5 .* K .* (cosh(2.*chi) + cos(2.*xi) - 2.*cosh(chi).*cos(xi));

    % A_n series (n = 1 to N): uses (n+1)² coefficient and cos(ξ)*cosh term
    % (instead of cosh(χ)*cosh as in stress_cc), reflecting the different
    % second-order derivatives of the stress function.
    for i = 1:N
        sig_xx = sig_xx + A(i) .* ( ...
            (i+1)^2  .* cosh((i+1).*chi) .* (cosh(chi) - cos(xi)) .* cos(i.*xi) ...
            -(i+1)   .* sinh((i+1).*chi) .* sinh(chi)              .* cos(i.*xi) ...
            + i      .* sin(xi)          .* sin(i.*xi)              .* cosh((i+1).*chi) ...
            + cos(i.*xi) .* cos(xi)      .* cosh((i+1).*chi) );
    end

    % B_1 constant term (same role as in stress_cc)
    sig_xx = sig_xx + B(1);

    % B_n series (n = 2 to N): analogous modification using (n-1)² and cos(ξ)
    for i = 2:N
        sig_xx = sig_xx + B(i) .* ( ...
            (i-1)^2  .* cosh((i-1).*chi) .* (cosh(chi) - cos(xi)) .* cos(i.*xi) ...
            -(i-1)   .* sinh((i-1).*chi) .* sinh(chi)              .* cos(i.*xi) ...
            + i      .* sin(xi)          .* sin(i.*xi)              .* cosh((i-1).*chi) ...
            + cos(i.*xi) .* cos(xi)      .* cosh((i-1).*chi) );
    end

    % Normalize by the focal distance c
    sig_xx = sig_xx / c;

end
