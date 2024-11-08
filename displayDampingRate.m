clear
clc

%% Generate plot showing damping rate as a function of \alpha and K_p
% optimisation - z = [alpha, kd, kp]
Q = @(z) [z(3), z(1)*z(3);
        z(1)*z(3), 1+z(1)^2*z(3)];
R = @(z) [z(1), 0;
        0, z(2)];
damping = @(z) min(eig(Q(z)*R(z)*Q(z)))/max(eig(Q(z)));

% Define the range for the first and second elements of z
alpha_range = 0.01:0.001:0.5;  % Adjust the range and step as necessary
kp_range = 0.1:0.1:25; 
kd = 10;

% Pre-allocate a matrix to store the damping values
damping_values2 = zeros(length(alpha_range), length(kp_range));

% Loop over the range of z1 and z2 values
for i = 1:length(alpha_range)
    for j = 1:length(kp_range)
        z = [alpha_range(i); kd; kp_range(j)];  % Assuming z(3) is constant, e.g., z(3) = 1
        damping_values2(i, j) = damping(z);
    end
end

% Create a meshgrid for plotting
[Z1, Z2] = meshgrid(alpha_range, kp_range);

% Plot the surface
figure;
% Generate a height map and contour plot of the objective function
contourf(Z1, Z2, damping_values2.', 100, 'LineStyle', 'none')
hold on
contour(Z1, Z2, damping_values2.', 10, 'k')
% Add labels to plot
xlabel('$$\alpha$$');
ylabel('$$K_p$$');
title('Contour plot of $$\frac{\lambda_{min}(QR_{cl}Q)}{\lambda_{max}(Q)}$$ as a function of $$\alpha$$ and $$K_p$$');