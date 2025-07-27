function dXdt = dynSysFastTime(X, epsilon, b, p, s)
% Extract variables
    phi = X(1);
    psi0 = X(2);
    A = X(3);
    B0 = X(4);
    
     % Define the first equation for d(phi)/dt (fast variable equation)
    num_phi = - A^2 * epsilon^2 *(b + epsilon + 2 * b * B0 * (-p + epsilon)) + ...
    epsilon * (-1 + b * B0 * epsilon) * (B0 * epsilon^5 + s * phi) + ...
    A * (b * B0 * (-1 + 2 * B0 * p) * epsilon^4 + B0 * (-1 + b - 2 * b * B0) * epsilon^5 + ...
    b * s * phi + s * epsilon * phi + epsilon^3*(1 + 2 * b * B0 - 2 * B0 * p +s * psi0));

    denom =  -A * (b + epsilon + b * B0 * (p - epsilon) * epsilon + B0 * p * epsilon ^ 2) + ...
    epsilon * (-1 + b * B0 * epsilon) * (-1 + B0 * epsilon * (-p + epsilon));

    dphidt = num_phi / (denom);
    
    % Define the second equation for d(psi0)/dt (fast variable equation)
    num_psi0 = A * B0 * (2 * p + 2 * B0 * p^2 * epsilon - epsilon^2) + 2 * b * B0 * (1 + B0 * (p - ...
        epsilon) * epsilon) * (- 1 + A + B0 * epsilon^2) + B0 * epsilon^2 * (1 + epsilon + ...
        B0 * epsilon * (p + (-1 + p) * epsilon)) + B0 * s * phi + s * (-1 + B0 * epsilon * (-p + epsilon)) * psi0;

    dpsi0dt = epsilon * num_psi0 / denom;

    % Define the rest of the slow equations
    dAdt = epsilon * phi;
    dB0dt = epsilon * psi0;

    % Return as a column vector
    dXdt = [dphidt; dpsi0dt; dAdt; dB0dt];

end


