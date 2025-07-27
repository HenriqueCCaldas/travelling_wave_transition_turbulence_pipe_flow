function plotSlowFlowVectorField(epsilon, b, p, s, ax)
    % Define a grid in Psi0, A, and B0 space. If needed, change the last
    % value to refine the mesh for the vector field
    [psi0, A, B0] = meshgrid(linspace(-0.2, 0.4, 16), ...
                             linspace(0, 1, 40), ...
                             linspace(0, 10, 16) ...
                             );
    
    psi0Flat = psi0(:);
    AFlat = A(:);
    B0Flat = B0(:);
    
    % Allocate derivative arrays
    dpsi0 = zeros(size(psi0Flat));
    dA = zeros(size(AFlat));
    dB0 = zeros(size(B0Flat));
    
    % Calculate derivatives at each point
    for i = 1:length(psi0Flat)
        Y = [psi0Flat(i); AFlat(i); B0Flat(i)];
        dY = slowFlow(Y, epsilon, b, p, s);
        dpsi0(i) = dY(1);
        dA(i) = dY(2);
        dB0(i) = dY(3);
    end
    
    % Compute the norm of each vector for scaling
    normFactor = sqrt(dpsi0.^2 + dA.^2 + dB0.^2);
    
    % Avoid division by zero or NaNs
    normFactor(normFactor < 1e-10 | isnan(normFactor)) = 1;
    
    % Scale vectors to uniform length for visualization
    scaleFactor = 0.05; % Adjust this value to control arrow length
    dpsi0 = scaleFactor * dpsi0 ./ normFactor;
    dA = scaleFactor * dA ./ normFactor;
    dB0 = scaleFactor * dB0 ./ normFactor;

    % create the color map and define the angls
    theta = atan2(dA, dpsi0);
    cmap = jet(256);
    
    if max(theta) == min(theta)
        % All angles are the same, use middle of colormap
        colorIndices = ones(size(theta)) * round(size(cmap, 1)/2);
    else
        % Map theta to indices in the colormap
        thetaNorm = (theta - min(theta)) / (max(theta) - min(theta));
        colorIndices = round(thetaNorm * (size(cmap, 1) - 1)) + 1;
        
        % Ensure indices are within valid range
        colorIndices = max(1, min(size(cmap, 1), colorIndices));
    end

    %verify is there is already an axis
    if nargin < 5 || isempty(ax)
        figure;
        ax = axes;
    end
    
    axes(ax);
    hold(ax, 'on');
    
    % Plot vectors as arrows
    for i = 1:length(psi0Flat)
        x = psi0Flat(i);
        y = AFlat(i);
        z = B0Flat(i);
        dx = dpsi0(i);
        dy = dA(i);
        dz = dB0(i);
        color = cmap(colorIndices(i), :);
        quiver3(x, y, z, dx, dy, dz, 0, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    end
    
    % Set labels and title
    xlabel('\psi_0', 'FontSize', 12);
    ylabel('A', 'FontSize', 12);
    zlabel('B_0', 'FontSize', 12);
    title('Slow Flow Vector Field', 'FontSize', 14);
    grid(ax, 'on');
    view(ax, 45, 30);
    colormap(ax, cmap);
    c = colorbar(ax);
    c.Label.String = 'Direction Angle';
end