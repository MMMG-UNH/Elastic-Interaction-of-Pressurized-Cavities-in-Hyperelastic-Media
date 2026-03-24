% =========================================================================
% stress_cx.m
%
% Computes the shear stress component σ_χξ in bipolar coordinates (χ, ξ).
%
% This is the off-diagonal stress acting on surfaces of constant χ and ξ.
% On the hole boundary (χ = χ₀), σ_χξ must be zero (traction-free in the
% tangential direction for a pressure-only loading).
%
% The expression implements Eq. (A8) from Saeedi & Kothari (2025), which
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
%   sig_cx - σ_χξ shear stress component at the given (chi, xi) coordinates
%            [same units as applied pressure P]
%
% Note: The shear stress formula is structurally simpler than σ_χχ and
% σ_ξξ. The leading K term involves sinh(χ)·sin(ξ), and all series terms
% involve only sin(nξ) (not cos(nξ)), reflecting the antisymmetry of shear
% about the symmetry planes (ξ = 0 and ξ = π).
% On both symmetry planes, sin(nξ) = 0, so σ_χξ = 0 there, as expected.
% =========================================================================

function sig_cx = stress_cx(chi, xi, A, B, K, N, c)

    % Leading term: −K·sinh(χ)·sin(ξ)
    % This term is zero at ξ=0 and ξ=π (symmetry planes) and antisymmetric
    % about these planes, consistent with the loading symmetry.
    sig_cx = -K .* (sinh(chi) .* sin(xi));

    % A_n series (n = 1 to N): each term has the factor n*(n+1) and
    % uses sinh((n+1)*chi) — note sinh not cosh, unlike σ_χχ and σ_ξξ.
    % The (cosh(chi) - cos(xi)) factor is the unnormalized Jacobian.
    for i = 1:N
        sig_cx = sig_cx + A(i) .* ( ...
            i .* (i+1) .* sinh((i+1).*chi) .* (cosh(chi) - cos(xi)) .* sin(i.*xi) );
    end

    % B_n series (n = 2 to N): analogous, using n*(n-1) and sinh((n-1)*chi).
    % No B_1 term here because the B_1 contribution to shear is zero
    % (the n=1 shear term vanishes: 1*(1-1) = 0).
    for i = 2:N
        sig_cx = sig_cx + B(i) .* ( ...
            i .* (i-1) .* sinh((i-1).*chi) .* (cosh(chi) - cos(xi)) .* sin(i.*xi) );
    end

    % Normalize by the focal distance c
    sig_cx = sig_cx / c;

end
