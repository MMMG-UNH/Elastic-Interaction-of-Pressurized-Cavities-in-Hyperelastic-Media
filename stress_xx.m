function sig_xx = stress_xx(eta, xi, A, B, K, N,a)
sig_xx = 0.5.*K.*(cosh(2.*eta)+cos(2.*xi) -2.*cosh(eta).*cos(xi));   
    for i=1:N
       sig_xx =  sig_xx + A(i).*((i+1)^2.*cosh((i+1).*eta).*(cosh(eta) - cos(xi)).*cos(i.*xi) -(i+1).*sinh((i+1).*eta).*sinh(eta).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i+1).*eta) + cos(i.*xi).*cos(xi).*cosh((i+1).*eta));
    end
    
    for i=2:N
        sig_xx = sig_xx +B(i).*((i-1)^2.*cosh((i-1).*eta).*(cosh(eta) - cos(xi)).*cos(i.*xi) -(i-1).*sinh((i-1).*eta).*sinh(eta).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i-1).*eta) + cos(i.*xi).*cos(xi).*cosh((i-1).*eta));
    end
    sig_xx = sig_xx + B(1);
    sig_xx = sig_xx/a;
end