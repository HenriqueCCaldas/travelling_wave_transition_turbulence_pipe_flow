function [t, X, Y] = heteroclinicOrbit(epsilon, b, p, s)

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

    % Initial time span for integration
    tspan = [0 3000];
    %Backward integration
    %tspan = [3000 0];
    % Initial guess for the trajectory

    % Perturb the initial trivial equilibrium (choose wheter you want
    % starting in trivial or nontrivial

    %X0 = X1 + 0.01*[0.091226;0;0.99583; 0];
    
    %X0 = X1 + 0.01 * [0.091226 ; 0; 0.99583 ; 0];

    X0 = X2 + 0.1 * [-0.00490214;0.0558524;-0.000371195;0.971917];
    Y0 = X0([2,3,4]);
    
    %X0 = X2 + 0.1*[-0.000895199 ; 0.824257 ; 0.0880286 ; 0.122488]; 
    % Gives good results for a heteroclinic conncetion between the 2 equilibria 
    
    % Solve the system with the initial condition
    [t, X] = ode15s(@(t, X) dynamicalSys(X, epsilon, b, p, s), tspan, X0);
    [~, Y] = ode15s(@(t, Y) slowFlow(Y,epsilon,b,p,s),tspan, Y0);
    
%{ 
Plot the trajectory in 2D (A vs B0)
figure;
%plot(X(:,3), X(:,4), 'b', 'LineWidth', 2);
xlim([-1 1])
ylim([0 10])
hold on;
plot(X1(3), X1(4), 'ro', 'MarkerSize', 10);
plot(X2(3), X2(4), 'go', 'MarkerSize', 10);
xlabel('A'); ylabel('B0');
title('Heteroclinic Orbit Exploration');
grid on;
traj = animatedline('Color', 'b', 'LineWidth', 2);
    for i = 1:length(t)
        addpoints(traj, X(i,3), X(i,4));
        drawnow;
        pause(0.03); 
    end
legend('Trivial Equilibrium', 'Nontrivial Equilibrium','Trajectory');
%}
    %Plot the trajectory in 3D
    figure;
    
    % Plot 3D trajectory
    plot3(X(:,2), X(:,3), X(:,4), 'b', 'LineWidth', 2);
    hold on;
    plot3(Y(:,1), Y(:,2), Y(:,3), 'r', 'LineWidth', 2);
    hold on;
    traj1 = animatedline('Color', 'b', 'LineWidth', 2);
    traj2 = animatedline('Color', 'r','LineWidth',2);

    % Plot equilibrium points in 3D
    plot3(X1(2), X1(3), X1(4), 'ro', 'MarkerSize', 10);
    plot3(X2(2), X2(3), X2(4), 'go', 'MarkerSize', 10);
    
    % Set axis labels and title
    xlim([-1 1])  % X-axis limits
    ylim([0 1])  % Y-axis limits
    zlim([0 10])  % Z-axis limits
    xlabel('Psi0'); 
    ylabel('A');
    zlabel('B0');
    title('Heteroclinic Orbit Exploration in 3D');

    % Adjust view and grid
    view(45, 30);  % Azimuth and elevation angles
    grid on;
    
    %{
    %Adding points to the trajetory
    for i = 1:length(t)
        %addpoints(traj1, X(i,2), X(i,3), X(i,4));
        addpoints(traj2, Y(i,1), Y(i,2), Y(i,3));
        drawnow;
        pause(0.03);
    end
    %}
    % Add legend
    legend('Trajectory', 'Trivial Equilibrium', 'Nontrivial Equilibrium');
end