function sing = singularity(i,X,b, epsilon, p)
    A = X(i,3);
    B0 = X(i,4);

    sing = - A * (b + epsilon + b * B0 * (p - epsilon) * epsilon + ...
                  B0 * p * epsilon^2) + ...
           epsilon * (-1 + b * B0 * epsilon) * ...
                   (-1 + B0 * epsilon * (-p + epsilon));

end