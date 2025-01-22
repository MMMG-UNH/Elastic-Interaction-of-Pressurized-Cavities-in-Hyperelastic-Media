function sig_ee = stress_ee(eta, xi, A, B, K, N,a)
sig_ee = -0.5.*K.*(cosh(2.*eta)+cos(2.*xi) -2.*cosh(eta).*cos(xi));   
    for i=1:N
       sig_ee =  sig_ee + A(i).*(-(i)^2.*cosh((i+1).*eta).*(cosh(eta) - cos(xi)).*cos(i.*xi) -(i+1).*sinh((i+1).*eta).*sinh(eta).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i+1).*eta) + cos(i.*xi).*cosh(eta).*cosh((i+1).*eta));
    end
    
    for i=2:N
        sig_ee = sig_ee +B(i).*(-(i)^2.*cosh((i-1).*eta).*(cosh(eta) - cos(xi)).*cos(i.*xi) -(i-1).*sinh((i-1).*eta).*sinh(eta).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i-1).*eta) + cos(i.*xi).*cosh(eta).*cosh((i-1).*eta));
    end
    sig_ee = sig_ee + B(1);
    sig_ee = sig_ee/a;
end