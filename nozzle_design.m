% Supersonic Wind Tunnel Nozzle Design (Method of Characteristics Geometry)
% Target Mach: ~1.64
% Method: Reflection/Cancellation (8 segments)
clear; clc; close all;

%% 1. INPUT PARAMETERS
gamma = 1.4;                % Specific heat ratio for Air
L_seg = 20;                 % Segment length in mm (2 cm)
theta_step = 2;             % Angle change per segment (degrees)
n_exp = 4;                  % Number of expansion segments
n_cancel = 4;               % Number of cancellation segments

% Define the Wall Angle Schedule (Theta)
% 4 segments ramping UP, 4 segments ramping DOWN to 0
theta_schedule = [2, 4, 6, 8, 6, 4, 2, 0]; 

%% 2. AERODYNAMIC CALCULATIONS
% The max wall angle (8 deg) represents half the total expansion in this method.
% Therefore, total Prandtl-Meyer angle (nu) at exit = 16 degrees.
theta_max = max(theta_schedule);
nu_exit_target = 2 * theta_max * (pi/180); % Convert to radians

% Solve for Design Mach Number using Inverse Prandtl-Meyer
% Function defined at bottom of script
M_exit = fzero(@(M) PM_Function(M, gamma) - nu_exit_target, [1.01, 3.0]);

% Calculate Area Ratio (A_exit / A_throat) based on Isentropic Flow
AR = (1/M_exit) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M_exit^2))^((gamma+1)/(2*(gamma-1)));

fprintf('------------------------------------------------\n');
fprintf('AERODYNAMIC DESIGN RESULTS\n');
fprintf('------------------------------------------------\n');
fprintf('Max Wall Angle:      %.2f deg\n', theta_max);
fprintf('Target P-M Angle:    %.2f deg\n', rad2deg(nu_exit_target));
fprintf('Calculated Mach:     %.4f\n', M_exit);
fprintf('Required Area Ratio: %.4f\n', AR);

%% 3. GEOMETRIC INTEGRATION
% Calculate total change in height (delta_y) from the wall segments
delta_y = 0;
for i = 1:length(theta_schedule)
    delta_y = delta_y + L_seg * sind(theta_schedule(i));
end

% Calculate necessary throat half-height
h_throat = delta_y / (AR - 1);
h_exit = h_throat * AR;

fprintf('Throat Half-Height:  %.2f mm\n', h_throat);
fprintf('Exit Half-Height:    %.2f mm\n', h_exit);
fprintf('Total Length:        %.2f mm\n', sum(L_seg * cosd(theta_schedule)));
fprintf('------------------------------------------------\n');

% Generate Coordinates
num_points = length(theta_schedule) + 1;
x = zeros(1, num_points);
y_upper = zeros(1, num_points);

% Starting Point (Throat)
x(1) = 0;
y_upper(1) = h_throat;

% Integrate downstream
for i = 1:length(theta_schedule)
    x(i+1) = x(i) + L_seg * cosd(theta_schedule(i));
    y_upper(i+1) = y_upper(i) + L_seg * sind(theta_schedule(i));
end

% Create Lower Wall (Symmetric)
y_lower = -y_upper;

%% 4. PLOTTING (CORRECTED)
figure('Name', 'Supersonic Nozzle Design', 'Color', 'w');
hold on; grid on; axis equal;

% Plot Centerline FIRST
plot([min(x)-10, max(x)+10], [0, 0], 'k-.', 'LineWidth', 1, 'DisplayName', 'Centerline');

% Plot Walls
plot(x, y_upper, 'b.-', 'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', 'Nozzle Walls');
plot(x, y_lower, 'b.-', 'LineWidth', 2, 'MarkerSize', 15, 'HandleVisibility', 'off'); % Hide from legend

% Draw Vertical Segment Lines (visual aid for MOC regions) - HIDDEN FROM LEGEND
for i = 1:length(x)
    line([x(i) x(i)], [y_lower(i) y_upper(i)], ...
        'Color', [0.8 0.8 0.8], ...
        'LineStyle', ':', ...
        'HandleVisibility', 'off'); % <--- This fixes the legend clutter
end

% Fill the shape for better visualization - HIDDEN FROM LEGEND
patch([x fliplr(x)], [y_upper fliplr(y_lower)], 'c', ...
    'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off'); % <--- This fixes the legend clutter

% Formatting
title(['MOC Nozzle Design (M = ' num2str(M_exit, '%.2f') ')']);
xlabel('Length (mm)');
ylabel('Height (mm)');
legend('Location', 'best'); % Now this will only show "Centerline" and "Nozzle Walls"
xlim([-10, max(x)+10]);
ylim([-h_exit*1.5, h_exit*1.5]);

%% 5. OUTPUT COORDINATES FOR CAD
fprintf('\nCOORDINATES FOR CAD (Upper Wall):\n');
fprintf('Point\t X (mm)\t\t Y (mm)\n');
fprintf('------------------------------\n');
for i = 1:length(x)
    fprintf('P%d\t %.3f\t\t %.3f\n', i-1, x(i), y_upper(i));
end

%% HELPER FUNCTION: PRANDTL-MEYER
function nu = PM_Function(M, g)
    % Calculates Prandtl-Meyer angle (in radians) for a given Mach Number
    term1 = sqrt((g+1)/(g-1));
    term2 = atan(sqrt((g-1)/(g+1) * (M^2 - 1)));
    term3 = atan(sqrt(M^2 - 1));
    nu = term1 * term2 - term3;
end