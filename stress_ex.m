function sig_ex = stress_ex(eta, xi, A, B, K, N,a)
sig_ex = -K.*(sinh(eta).*sin(xi));   
    for i=1:N
       sig_ex =  sig_ex + A(i).*(i.*(i+1).*sinh((i+1).*eta).*(cosh(eta) - cos(xi)).*sin(i.*xi));
    end
    
    for i=2:N
        sig_ex = sig_ex +B(i).*(i.*(i-1).*sinh((i-1).*eta).*(cosh(eta) - cos(xi)).*sin(i.*xi));
    end
    sig_ex = sig_ex/a;
end