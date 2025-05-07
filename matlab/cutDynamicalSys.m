function dXdt = cutDynamicalSys(X, epsilon, b, p, s)
% Extract variables
    phi = X(1);
    psi0 = X(2);
    A = X(3);
    B0 = X(4);
    
     % Define the first equation for d(phi)/dt (fast variable equation)
    dphidt = (-s * phi)/(epsilon) + (B0*p*s*phi) + epsilon * (A - 2*A*B0*p ...
        - B0*s*phi - B0^2 * p^2 * s * phi);
    
    % Define the second equation for d(psi0)/dt (slow variable equation)
    dpsi0dt = - (2*(-1+A)*b*B0+2*A*B0*p + B0 * s * phi -s * psi0)/(A*b) + ...
        epsilon*(2*b*B0-4*A*b*B0+2*A^2*b*B0-2*A*B0*p+2*A^2*B0*p-B0*s*phi + ...
        A * B0 *s * phi +A*b*B0^2*p*s*phi+s*psi0-A*s*psi0)/(A^2*b^2);

    % Define the rest of the slow equations
    dAdt = phi;
    dB0dt = psi0;

    % Return as a column vector
    dXdt = [dphidt; dpsi0dt; dAdt; dB0dt];

end