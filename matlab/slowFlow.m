function dXdt = slowFlow(Y, b, p, s)
%We specify the slow flow of the previous fast/slow system. 
    psi0 = Y(1);
    A = Y(2);
    B0 = Y(3);
    dpsi0 = -(2*(A-1)*b*B0+2*A*p*B0-s*psi0)/(A*b);
    dA = 0;
    dB0 = psi0;

    %Return as a collumn vector for the 3 dimensions

    dXdt = [dpsi0;dA;dB0];

end


