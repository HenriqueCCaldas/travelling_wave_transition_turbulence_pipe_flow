function dXdt = dynamicalSysOnlyS(X, epsilon, s)
% Extract variables
    phi = X(1);
    psi0 = X(2);
    A = X(3);
    B0 = X(4);
    
     % Define the first equation for d(phi)/dt (fast variable equation)
    num_phi = - (2* A^2 * (B0 - epsilon) * epsilon + A * epsilon * (epsilon + B0 * ...
        (-2 + 2 * B0 + (-2 + epsilon) * epsilon)) + (-1 + B0) * (B0 * epsilon^3 + ...
        s * phi) + A * s * (2 * phi + psi0));

    denom = epsilon * (-1 + B0^2 + A * (2 + 3 * B0));

    dphidt = num_phi / denom; 
    
    num_psi0 = - B0 * epsilon * (-2 + 2 * B0^2 + A * (6 + 10 * B0 - epsilon) + ...
        epsilon + epsilon^2 + B0 * epsilon * (1 + 2 * epsilon)) - B0 * s * phi + ...
        (1 + B0) * s * psi0;

    dpsi0dt = num_psi0 / denom;

    % Define the rest of the slow equations
    dAdt = phi;
    dB0dt = psi0;

    % Return as a column vector
    dXdt = [dphidt; dpsi0dt; dAdt; dB0dt];

end