function [t, X, Y] = heteroclinicOrbit(initalCond,onlyS,plotFullSystem,fullSystem,fastTimeScale,plotSlowFlow,fullSlowFlow,vecFieldSlowFlow, epsilon, b, p, s,tfinal)

    % Define initial conditions for the two equilibria
    % Trivial equilibrium
    X1 = [0; 0; 0; 0];

    if(onlyS == false)
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
    
    else

        Aeq = (4 - epsilon * (3 + 5 * epsilon) + ...
            sqrt(16 + (-1 + epsilon) * ...
            epsilon * (24 + epsilon * (7 +  epsilon))))/(24);

        Beq = (4 + (-1 + epsilon) * epsilon - sqrt(...
            16 + (-1 + epsilon) * epsilon * (24 + epsilon * (7 + ...
            epsilon))))/(8 * epsilon);
    end

    X2 = [0; 0; Aeq; Beq];

    %Nontrivial Vector

    
    %Nontrivial and nonphysical  equilibrium

    num3 = 4*b - 2*epsilon^2 - (2*b*epsilon^2)/p - 4*epsilon^3 - (2*b*epsilon^3)/p - ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+(-4*b*p-2*b*epsilon^2 + 2*p*epsilon^2 - 2*b*epsilon^3)^2)/p; 
    AeqV2 = num3 / (8*(b+p));

    num4 = 4*b*p + 2*b*epsilon^2 - 2*p*epsilon^2 + 2*b*epsilon^3 + ...
           sqrt(-16*b*p*epsilon^2*(2*b-epsilon^2-epsilon^3)+...
           (-4*b*p-2*b*epsilon^2+2*p*epsilon^2 - 2*b*epsilon^3)^2);

    BeqV2 = num4 / (8 *b*p*epsilon^2);

    % Initial time span for integration
    tspan = [0 tfinal];

    %Backward integration
    %tspan = [3000 0];
    % Perturb the initial trivial equilibrium (choose wheter you want
    % starting in trivial or nontrivial by choosing X0 or X1
    if(initalCond == "trivial")
        %Testing for shifthing from unstable direction of the trivial equilibria
        %for b = 0.04
        X0 = X1 + 0.01 * [0.0487508;0.0;0.998811;0.0];
        %Testing for shifting from stable direction of the trivial equilibria
        %for b = 0.04
        %X0 = X1 + 1 * [-0.898668;0;0.438629;0];
    else
        %Testing for shifting from unstable direction of the non-trivial
        %equilibrum

        if(b == 0.05)
            if(s == 0.08)
                X0 = X2 + 0.01 * [-0.007; 0.098; -0.001; 0.953];
            else
                %irrelevant scenario
                X0 = X2 + 0.01 * [-0.00277756;0.0126069;-0.000608936;0.991606];
            end
        elseif(b == 0.10)

            if(s == 0.08)
                X0 = X2 + 0.01 * [0.01; -0.143; -0.001;- 0.912];
            else
                X0 = X2 + 0.01 * [0.007; -0.042; 0.001;- 0.967];
            end
        else
            X0 = X2 + 0.01 * [0.004; -0.981; 0.001;- 0.189];
        end
    end

    Y0 = X0([2,3,4]);
    
    if(onlyS == true)
        [t,X] = ode15s(@(t,X) dynamicalSysOnlyS(X,epsilon,s), tspan, X0);
    else
        if(fullSystem == true)
            % Solve the system with the initial condition for different setups
            if(fastTimeScale == true)
                [t, X] = ode15s(@(t, X) dynSysFastTime(X, epsilon, b,p,s), tspan, X0);
            else
                %Option 2: Normal full dynamical system
                [t, X] = ode15s(@(t, X) dynamicalSys(X, epsilon, b, p, s), tspan, X0);
            end
        else
            %Option 2:Cut dynamical system up to epsilon^2
            [t, X] = ode15s(@(t, X) cutDynamicalSys(X, epsilon, b, p, s), tspan, X0);
        end
    end
    
    

    if(fullSlowFlow == true)
            %Option4: normal slow flow 
            [tY, Y] = ode15s(@(tY, Y) slowFlow(Y, epsilon, b, p, s), tspan, Y0);
    else
            %Option 3: quasi slow flow flow with eps terms
            [tY, Y] = ode15s(@(t, Y) quasiSlowFlow(Y, X0(1),t,epsilon,b,p,s),tspan,Y0);
    end
        
    % Sort time 
    [t, idx] = sort(t);
    X = X(idx, :);
    [tY, idy] = sort(tY);
    Y = Y(idy, :);
    
    fig = figure;
    mainLayout = tiledlayout(fig,4,7);

    % 3D Plot on the left tile
    ax = nexttile(mainLayout,[4,3]);
    hold (ax,'on');

    plot3(ax,X1(2), X1(3), X1(4), 'ro', 'MarkerSize', 10); % Trivial
    plot3(ax,X2(2), X2(3), X2(4), 'go', 'MarkerSize', 10); % Nontrivial
    %plot3(ax,X3(2), X3(3), X3(4), 'mo', 'MarkerSize', 10); % Nonphysical
    %limits of the plot
    
    xlim(ax,[-20, 20]);
    ylim(ax,[-1, 1]);
    zlim(ax,[-10 30]);
    xlabel(ax,'\psi_0','Interpreter','tex'); 
    ylabel("A",'Interpreter','tex'); 
    zlabel("B_0",'Interpreter','tex');
    if(fastTimeScale == true)
        title(ax,['Fast time scale - Phase space for \psi_{0} , A and B_0 with b = ' num2str(b) ' and s = ' num2str(s)],'Interpreter','tex');
    else
        title(ax,['Phase space for \psi_{0} , A and B_0 (b = ' num2str(b) ' and s = ' num2str(s) ')'],'Interpreter','tex');
    end
    view(ax, 45, 15); 
    grid (ax,'on');
    ax.FontSize = 20;

    % Add trajectory points
    N = min(length(t), size(X,1));
    M = min(length(tY), size(Y,1));

    %Slow flow vector field 
    if (plotSlowFlow == true)
        if (vecFieldSlowFlow == false)            
            traj2 = animatedline('Color', 'b', 'LineWidth', 2);
            hold on;
            for i = 1:M
            if any(isnan(Y(i,:))) || any(isinf(Y(i,:))), continue; end
            addpoints(traj2, Y(i,1), Y(i,2), Y(i,3));
            %fprintf('Y at t = %.5f: [%.4f, %.4f, %.4f]\n', t(i), Y(i,1), Y(i,2), Y(i,3));
            drawnow;
            pause(0.0001);
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
            singAtT = singularity(i,X, b, epsilon, p);
            fprintf('X at t = %.5f: [%.4f,%.4f, %.4f, %.4f]. Singularity: %.20f\n', t(i),X(i,1), X(i,2), X(i,3), X(i,4),singAtT);
            drawnow;
            pause(0.001);
        end
    end
    if (plotSlowFlow == true)
        if plotFullSystem == true
            legend(ax,'Trivial Eq.','Nontrivial Eq. 1','Slow flow','Full system');
        else
            legend(ax,'Trivial Eq.','Nontrivial Eq. 1','Slow flow');
        end
    else
        if plotFullSystem == true
            legend(ax,'Trivial Eq.','Nontrivial Eq. 1','Full system');
        else
            legend(ax,'Trivial Eq.','Nontrivial Eq. 1');
        end
    end

   % Subplot 1:Phi
   nexttile(mainLayout,[1,4])  
   plot(t, X(:,1));
   title('Time Series of \phi in the full system');
   ylabel('\phi');
    
   % Subplot 2: Psi_0
   nexttile(mainLayout,[1,4])
   plot(t, X(:,2));
   title('Time Series of \psi_0 in the full system');
   ylabel('\psi_0');
    
   % Subplot 3: A
   nexttile(mainLayout,[1,4])
   plot(t, X(:,3));
   title('Time Series of A in the full system');
   ylabel('A');
    
   % Subplot 4: B_0
   nexttile(mainLayout,[1,4])
   plot(t, X(:,4));
   title('Time Series of B_0 in the full system');
   xlabel('Time');
   ylabel('B_0');

   allAxes = findall(fig, 'Type', 'Axes');
   for ax = allAxes'
    ax.FontSize = 20;
   end

end