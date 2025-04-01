function dXdt = slowFlow(Y, epsilon, b, p, s)
%We specify the slow flow of the previous fast/slow system. For this
%representation, we consider the terms up until order epsilon, to correct
%for the higher order terms.
    psi0 = Y(1);
    A = Y(2);
    B0 = Y(3);
    dpsi0 = -(2*(A-1)*b*B0+2*A*p*B0-s*psi0)/(A*b) + ...
            epsilon*(2*b*B0*A^2+2*B0*p*A^2-4*b*B0*A-2*B0*p*A-s*psi0*A+2*b*B0+s*psi0)/(A^2*b^2);
    %dpsi0 = -(2*(A-1)*b*B0+2*A*p*B0-s*psi0)/(A*b);
    dA = 0;
    dB0 = psi0;

    %Return as a collumn vector for the 3 dimensions

    dXdt = [dpsi0;dA;dB0];

end

