function sig_cx = stress_cx(chi, xi, A, B, K, N, c)
sig_cx = -K.*(sinh(chi).*sin(xi));   
    for i=1:N
       sig_cx =  sig_cx + A(i).*(i.*(i+1).*sinh((i+1).*chi).*(cosh(chi) - cos(xi)).*sin(i.*xi));
    end
    
    for i=2:N
        sig_cx = sig_cx +B(i).*(i.*(i-1).*sinh((i-1).*chi).*(cosh(chi) - cos(xi)).*sin(i.*xi));
    end
    sig_cx = sig_cx/c;
end