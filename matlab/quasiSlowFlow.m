function dXdt = quasiSlowFlow(Y, phi0, t, epsilon, b, p, s)
%We specify the slow flow of the previous fast/slow system. For this
%representation, we consider the terms up until order epsilon, to correct
%for the higher order terms.
    %Assign the current values of the relevant slow variables
    psi0 = Y(1);
    A = Y(2);
    B0 = Y(3);
    %Assign the current evaluation of psi (phi0 is obtained as an input)
    phi = phi0 * exp( (- s/epsilon + B0 * p * s) * t);

    dpsi0 = -(2 * ( A - 1) * b * B0 + 2 * A * p * B0 + B0 * s * phi- s * psi0)/(A*b) + ...
            epsilon * (2 * b * B0 - 4 * A * b * B0 + 2 * A^2 * b * B0 - 2 * A * B0 * p + 2 * A^2 * B0 * p - ...
            B0 * s * phi + A * b * B0^2 * p * s * phi + s * psi0 - A * s * psi0)/(A^2 * b^2);
    %dpsi0 = -(2*(A-1)*b*B0+2*A*p*B0-s*psi0)/(A*b);
    dA = phi;
    dB0 = psi0;

    %Return as a collumn vector for the 3 dimensions

    dXdt = [dpsi0;dA;dB0];

end


