function [t, X, Y] = heteroclinicOrbit(initalCond,plotFullSystem, fullSystem,plotSlowFlow,fullSlowFlow,vecFieldSlowFlow, epsilon, b, p, s)

    % Define initial conditions for the two equilibria
    % Trivial equilibrium
    X1 = [0; 0; 0; 0];
    
    % Definition of A in nontrivial equilibrium

    num1 = 4*b - 2*epsilon^2 - (2*b*epsilon^2)/p - 4*epsilon^3 - (2*b*epsilon^3)/p + ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+(-4*b*p-2*b*epsilon^2 + 2*p*epsilon^2 - 2*b*epsilon^3)^2)/p;
    den1 = 8*(b+p);
    Aeq = num1/den1;
    
    % Definition of B in nontrivial equilibrium

    num2 = 4*b*p + 2*b*epsilon^2 - 2*p*epsilon^2 + 2*b*epsilon^3 - ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+...
           (-4*b*p-2*b*epsilon^2+2*p*epsilon^2 - 2*b*epsilon^3)^2);
    den2 = 8 *b*p*epsilon^2;
    Beq = num2/(den2);

    %Nontrivial Vector

    X2 = [0; 0; Aeq; Beq];


    %Nontrivial and nonphysical  equilibrium

    num3 = 4*b - 2*epsilon^2 - (2*b*epsilon^2)/p - 4*epsilon^3 - (2*b*epsilon^3)/p - ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+(-4*b*p-2*b*epsilon^2 + 2*p*epsilon^2 - 2*b*epsilon^3)^2)/p; 
    AeqV2 = num3/den1;

    num4 = 4*b*p + 2*b*epsilon^2 - 2*p*epsilon^2 + 2*b*epsilon^3 + ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+...
           (-4*b*p-2*b*epsilon^2+2*p*epsilon^2 - 2*b*epsilon^3)^2);

    BeqV2 = num4/den2;

    X3 = [0; 0; AeqV2; BeqV2];

    % Initial time span for integration
    tspan = [0 3000];

    %Backward integration
    %tspan = [3000 0];
    % Perturb the initial trivial equilibrium (choose wheter you want
    % starting in trivial or nontrivial by choosing X0 or X1
    if(initalCond == "trivial")
        %Testing for shifthing from unstable direction of the trivial equilibria
        %for b = 0.04
        X0 = X1 + 1 * [0.0487508;0;0.998811;0];
        %Testing for shifting from stable direction of the trivial equilibria
        %for b = 0.04
        %X0 = X1 + 1 * [-0.898668;0;0.438629;0];
    else
        %Testing for shifting from unstable direction of the non-trivial equilibria
        %for b = 0.04 this yields a chaotic atractor of the system. This
        %presents an interesting chaotic effect. Most interesting effect for
        %now, as well as to show that the whole system some how shows the whole
        %dynamics
        X0 = X2 + 0.01 * [-0.00277756;0.0126069;-0.000608936;0.991606];
    
        %X0 =  X2 + 0.1*[-0.00277756; 0.0126069; -0.000608936 ; 0.991606];
    end

    Y0 = X0([2,3,4]);
    if(fullSystem == true)
        % Solve the system with the initial condition for different setups
        %Option 1:Normal full dynamical system
        [t, X] = ode15s(@(t, X) dynamicalSys(X, epsilon, b, p, s), tspan, X0);
    else
        %Option 2:Cut dynamical system up to epsilon^2
        [t, X] = ode15s(@(t, X) cutDynamicalSys(X, epsilon, b, p, s), tspan, X0);
    end

    if(fullSlowFlow == true)
        %Option4: normal slow flow 
        [tY, Y] = ode15s(@(tY, Y) slowFlow(Y, epsilon, b, p, s), tspan, Y0);
    else
        %Option 3: quasi slow flow flow with eps terms
        [tY, Y] = ode15s(@(t, Y) quasiSlowFlow(Y, X0(1),t,epsilon,b,p,s),tspan,Y0);
    end

    % Sort time (fixes weird jumps)
    [t, idx] = sort(t);
    X = X(idx, :);
    [tY, idy] = sort(tY);
    Y = Y(idy, :);

    % 3D Plot
    fig = figure;
    ax = axes('Parent', fig);
    hold on;

    plot3(ax,X1(2), X1(3), X1(4), 'ro', 'MarkerSize', 10); % Trivial
    plot3(ax,X2(2), X2(3), X2(4), 'go', 'MarkerSize', 10); % Nontrivial
    %plot3(ax,X3(2), X3(3), X3(4), 'mo', 'MarkerSize', 10); % Nonphysical
    %limits of the plot
    
    xlim([-20, 20]);
    ylim([-0.1, 1]);
    zlim([0 40]);
    xlabel('\phi_0','Interpreter','tex'); ylabel("A",'Interpreter','tex'); zlabel("B_0",'Interpreter','tex');
    title('Phase Space for \psi_{0} , A and B_0 ','Interpreter','tex');
    view(45, 30); grid on;

    % Add trajectory points
    N = min(length(t), size(X,1));
    M = min(length(tY), size(Y,1));

    %Uncomment for Slow flow vector field 
    if (plotSlowFlow == true)
        if (vecFieldSlowFlow == false)            
            traj2 = animatedline('Color', 'b', 'LineWidth', 2);
            hold on;
            for i = 1:M
            if any(isnan(Y(i,:))) || any(isinf(Y(i,:))), continue; end
            addpoints(traj2, Y(i,1), Y(i,2), Y(i,3));
            fprintf('Y at t = %.5f: [%.4f, %.4f, %.4f]\n', t(i), Y(i,1), Y(i,2), Y(i,3));
            drawnow;
            pause(0.001);
            end
        else
            plotSlowFlowVectorField(epsilon,b,p,s,ax);
        end
    end
    if (plotFullSystem == true)
        traj1 = animatedline(ax,'Color', 'r', 'LineWidth', 2);
        hold on;
        for i = 1:N
            if any(isnan(X(i,:))) || any(isinf(X(i,:))), continue; end
            addpoints(traj1, X(i,2), X(i,3), X(i,4));
            fprintf('X at t = %.5f: [%.4f, %.4f, %.4f]\n', t(i), X(i,2), X(i,3), X(i,4));
            drawnow;
            pause(0.0001);
        end
    end
    if (plotSlowFlow == true)
        if plotFullSystem == true
            legend('Trivial Eq.','Nontrivial Eq. 1','Slow Flow','Trajectory');
        else
            legend('Trivial Eq.','Nontrivial Eq. 1','Slow Flow');
        end
    else
        if plotFullSystem == true
            legend('Trivial Eq.','Nontrivial Eq. 1','Trajectory');
        else
            legend('Trivial Eq.','Nontrivial Eq. 1');
        end
    end

end