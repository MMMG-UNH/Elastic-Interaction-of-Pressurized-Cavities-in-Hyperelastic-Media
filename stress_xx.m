function sig_xx = stress_xx(chi, xi, A, B, K, N, c)
sig_xx = 0.5.*K.*(cosh(2.*chi)+cos(2.*xi) -2.*cosh(chi).*cos(xi));   
    for i=1:N
       sig_xx =  sig_xx + A(i).*((i+1)^2.*cosh((i+1).*chi).*(cosh(chi) - cos(xi)).*cos(i.*xi) -(i+1).*sinh((i+1).*chi).*sinh(chi).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i+1).*chi) + cos(i.*xi).*cos(xi).*cosh((i+1).*chi));
    end
    
    for i=2:N
        sig_xx = sig_xx +B(i).*((i-1)^2.*cosh((i-1).*chi).*(cosh(chi) - cos(xi)).*cos(i.*xi) -(i-1).*sinh((i-1).*chi).*sinh(chi).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i-1).*chi) + cos(i.*xi).*cos(xi).*cosh((i-1).*chi));
    end
    sig_xx = sig_xx + B(1);
    sig_xx = sig_xx/c;
end