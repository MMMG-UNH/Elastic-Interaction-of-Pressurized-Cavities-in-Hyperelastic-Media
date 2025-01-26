function sig_cc = stress_cc(chi, xi, A, B, K, N, c)
sig_cc = -0.5.*K.*(cosh(2.*chi)+cos(2.*xi) -2.*cosh(chi).*cos(xi));   
    for i=1:N
       sig_cc =  sig_cc + A(i).*(-(i)^2.*cosh((i+1).*chi).*(cosh(chi) - cos(xi)).*cos(i.*xi) -(i+1).*sinh((i+1).*chi).*sinh(chi).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i+1).*chi) + cos(i.*xi).*cosh(chi).*cosh((i+1).*chi));
    end
    
    for i=2:N
        sig_cc = sig_cc +B(i).*(-(i)^2.*cosh((i-1).*chi).*(cosh(chi) - cos(xi)).*cos(i.*xi) -(i-1).*sinh((i-1).*chi).*sinh(chi).*cos(i.*xi) ...
           +i.*sin(xi).*sin(i.*xi).*cosh((i-1).*chi) + cos(i.*xi).*cosh(chi).*cosh((i-1).*chi));
    end
    sig_cc = sig_cc + B(1);
    sig_cc = sig_cc/c;
end