function dXdt = dynamicalSys(X, epsilon, b, p, s)
% Extract variables
    phi = X(1);
    psi0 = X(2);
    A = X(3);
    B0 = X(4);
    
     % Define the first equation for d(phi)/dt (fast variable equation)
    num_phi = -A^2 * epsilon^2 * (2 * b * B0 * (epsilon - p) + b + epsilon) + ...
              A * (b * B0 * epsilon^4 * (2 * B0 * p - 1) + epsilon^3 * (2 * b * B0 - 2 * B0 * p + 1) + ...
              B0 * epsilon^5 * (-2 * b * B0 + b - 1) + b * s * phi + s * psi0 * epsilon^2 + s * epsilon * phi) + ...
              epsilon * (b * B0 * epsilon - 1) * (B0 * epsilon^5 + s * phi);
    
    denom = epsilon * (b * B0 * epsilon - 1) * (B0 * epsilon * (epsilon - p) - 1) - ...
            A * (b * B0 * epsilon * (p - epsilon) + b + B0 * p * epsilon^2 + epsilon);
    
    dphidt = num_phi / (denom * epsilon);
    
    % Define the second equation for d(psi0)/dt (fast variable equation)
    num_psi0 = 2 * b * B0 * epsilon * (A + B0 * epsilon^2 - 1) * (B0 * epsilon * (p - epsilon) + 1) + ...
               A * B0 * epsilon * (2 * B0 * p^2 * epsilon + 2 * p - epsilon^2) + ...
               s * psi0 * (B0 * epsilon * (epsilon - p) - 1) + ...
               B0 * epsilon^3 * (B0 * epsilon * ((p - 1) * epsilon + p) + epsilon + 1) + ...
               B0 * s * epsilon * phi;
    dpsi0dt = num_psi0 / denom;

    % Define the rest of the slow equations
    dAdt = phi;
    dB0dt = psi0;

    % Return as a column vector
    dXdt = [dphidt; dpsi0dt; dAdt; dB0dt];

end