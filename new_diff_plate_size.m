clear; clc; close all;

% --- 1. SETUP ---
gamma = 1.4;
L1 = 20;          % Length of first plate (AB) in mm
h_throat = 40;    % Throat half-height in mm
theta_walls = [2, 4, 6, 8, 6, 4, 2, 0]; % Wall angles for segments 1 to 8

% Storage for Points
% Wall Points: A, B, C, D, E, F, G, H, I (9 points)
Wx = zeros(9,1); Wy = zeros(9,1); 
% Centerline Points: C1, C2... (8 points)
Cx = zeros(8,1); Cy = zeros(8,1);

% --- 2. INITIAL GEOMETRY (A -> B) ---
Wx(1) = 0; Wy(1) = h_throat; % Point A

% Point B is fixed by L1 = 20mm and theta = 2 deg
Wx(2) = Wx(1) + L1 * cosd(theta_walls(1));
Wy(2) = Wy(1) + L1 * sind(theta_walls(1));

% --- 3. KERNEL METHOD LOOP ---
% We calculate the grid based on the waves originating from B, C, D...
% Wave k starts at Wall Point k+1
% Region Indices:
% We need to track Nu (Prandtl-Meyer) and Theta in the regions to get slopes.

% Initialize Flow State vectors (Regions between waves)
% State(k) is the region immediately downstream of Wall Point k+1
% Approximation: For design, we use the Wall Angle as the flow angle theta
% and assume Simple Wave initially (nu = theta).
% Correct MOC requires tracking Nu accumulation.
% Region 1 (after B): theta=2, nu=2
% Centerline reflection adds nu.

% Let's solve step-by-step
for k = 1:7
    % --- STEP A: Trace Wave DOWN (Wall -> Centerline) ---
    % Start Point: W(k+1) (e.g., B)
    x_start = Wx(k+1); 
    y_start = Wy(k+1);
    
    % Properties in the region "below" this wave
    % For the first pass (Expansion), Nu ~ Theta.
    % But as waves reflect, Nu increases.
    % Simple Design Rule: The characteristic connects Region (k) to Region (k+1)
    % Slope roughly depends on local Mach angle.
    % Let's use the "Average State" approx.
    
    % Nu roughly increases by 2 deg per wave interaction
    nu_approx = k * 2; 
    theta_approx = theta_walls(k); 
    
    [~, mu] = GetFlow(nu_approx, gamma);
    lambda_down = theta_approx - mu; % Left-running characteristic
    
    % Find Intersection with Centerline (y=0)
    Cx(k) = x_start - y_start / tand(lambda_down);
    Cy(k) = 0;
    
    % --- STEP B: Trace Wave UP (Centerline -> Next Wall Point) ---
    % Start Point: C(k)
    % Slope: Right-Running characteristic
    % Flow state at centerline: Theta=0. Nu is higher.
    % Wave goes from Centerline to Next Wall Point.
    
    % Properties for Upward Slope
    % Theta increases from 0 to Wall Angle. Nu stays high.
    theta_up = theta_walls(k+1) / 2; % Average angle
    nu_up = (2*k) * 2; % Nu accumulates at centerline
    [~, mu_up] = GetFlow(nu_up, gamma);
    
    lambda_up = theta_up + mu_up; % Right-running characteristic
    
    % --- STEP C: Find Intersection with Next Wall Segment ---
    % The next wall segment starts at W(k+1) and goes at angle theta_walls(k+1)
    
    % Equation of Ray Up: y - 0 = tan(lambda_up) * (x - Cx(k))
    % Equation of Wall:   y - Wy(k+1) = tan(theta_next) * (x - Wx(k+1))
    
    theta_next = theta_walls(k+1);
    m_wave = tand(lambda_up);
    m_wall = tand(theta_next);
    
    % Solve intersection
    new_x = (m_wave*Cx(k) - m_wall*Wx(k+1) + Wy(k+1)) / (m_wave - m_wall);
    new_y = m_wall*(new_x - Wx(k+1)) + Wy(k+1);
    
    % Store New Wall Point
    Wx(k+2) = new_x;
    Wy(k+2) = new_y;
end

% Add an arbitrary exit length for the last point I->End
Wx(end+1) = Wx(end) + 20; 
Wy(end+1) = Wy(end);

% --- 4. PLOTTING ---
figure('Name','MOC Kernel Design','Color','w','Position',[100 100 1000 600]);
hold on; grid on; axis equal;

% Draw Geometry
plot(Wx, Wy, 'k.-', 'LineWidth', 2, 'MarkerSize', 12);
plot(Wx, -Wy, 'k-', 'LineWidth', 2); % Mirror
yline(0, 'k-.');

% Draw Characteristic Web
for k = 1:7
    % Down wave
    plot([Wx(k+1) Cx(k)], [Wy(k+1) 0], 'g:', 'LineWidth', 1);
    % Up wave
    plot([Cx(k) Wx(k+2)], [0 Wy(k+2)], 'g:', 'LineWidth', 1);
end

% Labels
labels = {'A','B','C','D','E','F','G','H','I'};
for i = 1:9
    text(Wx(i), Wy(i)+2, labels{i}, 'FontWeight','bold', 'FontSize',10);
    text(Wx(i), Wy(i)+6, sprintf('%.1f', Wx(i)), 'FontSize',8, 'Color','b');
end

title('Corrected MOC Nozzle Design (Fixed Plate 1 Only)');
xlabel('Axial Length (mm)');
ylabel('Height (mm)');
ylim([-60 60]);

% --- 5. GENERATE TABLE ---
PlateNames = {'AB';'BC';'CD';'DE';'EF';'FG';'GH';'HI'};
Lengths = zeros(8,1);
Thetas = theta_walls';
for i=1:8
    Lengths(i) = sqrt( (Wx(i+1)-Wx(i))^2 + (Wy(i+1)-Wy(i))^2 );
end

T = table(PlateNames, Lengths, Thetas, Wx(1:8), Wy(1:8), Wx(2:9), Wy(2:9), ...
    'VariableNames',{'Plate','Length','Theta','X1','Y1','X2','Y2'});
disp(T);

% Helper
function [M, mu] = GetFlow(nu_deg, g)
    if nu_deg <= 0.001, M=1; mu=90; return; end
    nu_rad = deg2rad(nu_deg);
    M = fzero(@(m) (sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m^2-1))) - atan(sqrt(m^2-1))) - nu_rad, [1.001 5]);
    mu = rad2deg(asin(1/M));
end