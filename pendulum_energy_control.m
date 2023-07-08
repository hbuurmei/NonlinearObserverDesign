% Energy control method applied to pendulum
% with pendulum written as a control affine system
% NOTE: theta = 0 is defined to be verticaly upwards here

% Define system parameters
sys.b = 0.1; % damping
sys.m = 1; % mass
sys.g = 9.8; % gravity
sys.l = 1; % length
sys.u_max = 0.05 * sys.m*sys.g*sys.l; % max torque input

% Define sim parameters
eps = 1e-2; % tolerance
dt = 0.01;
t = 0;
max_time = 50;

% Simulation loop
x = [5*pi/4; 0];
x_des = zeros(2,1);
E0 = 0; % desired energy level (defined to be zero)
x_hist = [];
u_hist = [];
t_hist = [];

figure
hold on
grid on
h_origin = plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % filled point at origin
h_pendulum = plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % plotting the point mass
line_handle = line([0, 0], [0, -sys.l], 'LineWidth', 2); % plotting the pendulum
trajectory_line = plot(-sys.l * sin(x(1)), sys.l * cos(x(1)), '--k'); % dashed trajectory line
xlim([-sys.l-0.2, sys.l+0.2]); % Set fixed limits for the x-axis
ylim([-sys.l-0.2, sys.l+0.2]); % Set fixed limits for the y-axis
axis equal; % Set equal scaling for x and y axes
xlabel('x');
ylabel('y');

while norm(x-x_des) > eps
    E = system_energy(x, sys);
    u = -sys.u_max*extended_sign((E - E0)*x(2));
    x_dot = f(x, u, sys);
    x = x + dt * x_dot;
    x(1) = wrapTo2Pi(x(1));
    % fprintf('theta (deg): %i\n', x(1)*180/pi)
    % fprintf('u: %i\n', u)
    % fprintf('t: %i\n', t)
    t = t + dt;
    x_hist = [x_hist, x];
    u_hist = [u_hist, u];
    t_hist = [t_hist, t];
    if t > max_time
        break
    end

    % Update the plot
    x_pendulum = -sys.l * sin(x(1));
    y_pendulum = sys.l * cos(x(1));
    set(h_origin, 'XData', 0, 'YData', 0);
    set(h_pendulum, 'XData', x_pendulum, 'YData', y_pendulum);
    set(line_handle, 'XData', [0, x_pendulum], 'YData', [0, y_pendulum]);

    % Update trajectory line
    set(trajectory_line, 'XData', -sys.l*sin(x_hist(1, :)), 'YData', sys.l * cos(x_hist(1, :)));

    % Restore initial x-axis limits
    xlim([-sys.l-0.2, sys.l+0.2]);
    drawnow;
end

hold off

%% Plot results
figure
plot(t_hist, x_hist(1, :))
figure
plot(t_hist, u_hist)

%% Helper functions
function E = system_energy(x, sys)
    % Energy of system
    theta = x(1);
    theta_dot = x(2);
    E = 0.5*sys.m*sys.l^2*theta_dot^2 + sys.m*sys.g*sys.l*(cos(theta)-1);
end

function x_dot = f(x, u, sys)
    % Dynamics   
    theta = x(1);
    theta_dot = x(2);
    theta_ddot = 1/(sys.m*sys.l^2) * (-sys.b*theta_dot + sys.m*sys.g*sys.l*sin(theta) + u);
    x_dot = [theta_dot; theta_ddot];
end

function sgn = extended_sign(z)
    % Sign function but equals 1 if input is zero
    if z >= 0
        sgn = 1;
    else
        sgn = -1;
    end
end
