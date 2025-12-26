clear; clc; close all;

% --- 1. INPUT PARAMETERS ---
gamma = 1.4;
L_fixed = 20;   % mm (Fixed length for Expansion Plates)
h_throat = 40;  % mm (Throat Half-Height)
theta_step = 2; % deg
n_waves = 4;    % Number of waves

% --- 2. GENERATE EXPANSION WALL (A -> B -> C -> D -> E) ---
% We fix the first 4 segments to 20mm to create the "Kernel"
xw = zeros(1, 2*n_waves + 1);
yw = zeros(1, 2*n_waves + 1);

% Point A
xw(1) = 0; yw(1) = h_throat;

% Calculate B, C, D, E
current_theta = 0;
for i = 1:n_waves
    current_theta = current_theta + theta_step;
    xw(i+1) = xw(i) + L_fixed * cosd(current_theta);
    yw(i+1) = yw(i) + L_fixed * sind(current_theta);
end

% --- 3. MOC GRID SOLVER ---
Px = zeros(n_waves, n_waves);
Py = zeros(n_waves, n_waves);
Cx = zeros(1, n_waves); % Centerline

% We solve Reflected Wave by Reflected Wave (Column j)
for j = 1:n_waves
    
    % --- A. Find Centerline Point C(j) ---
    % Intersection of Left-Running Wave with y=0
    if j == 1
        src_x = xw(2); src_y = yw(2); % Start at B
        theta_loc = 1 * theta_step;
        nu_loc = theta_loc;
    else
        src_x = Px(j, j-1); src_y = Py(j, j-1); % Start at P(j, j-1)
        theta_loc = theta_step; 
        nu_loc = (2*j - 1) * theta_step;
    end
    
    [~, mu] = GetFlow(nu_loc, gamma);
    lambda = theta_loc - mu;
    Cx(j) = src_x - src_y / tand(lambda);
    
    % --- B. Internal Grid Points P(i,j) ---
    for i = j:n_waves
        % Upper Source (Left Running)
        if j == 1
            xu = xw(i+1); yu = yw(i+1); % From Wall
        else
            xu = Px(i, j-1); yu = Py(i, j-1);
        end
        
        % Lower Source (Right Running)
        if i == j
            xl = Cx(j); yl = 0; % From Centerline
        else
            xl = Px(i-1, j); yl = Py(i-1, j);
        end
        
        % Region Properties (Target Cell)
        theta_cell = (i - j) * theta_step;
        nu_cell = (i + j) * theta_step;
        [~, mu_c] = GetFlow(nu_cell, gamma);
        
        % Slopes
        m_minus = tand(theta_cell - mu_c);
        m_plus  = tand(theta_cell + mu_c);
        
        % Intersect
        x_int = (m_minus*xu - m_plus*xl + yl - yu) / (m_minus - m_plus);
        y_int = m_minus*(x_int - xu) + yu;
        
        Px(i,j) = x_int;
        Py(i,j) = y_int;
    end
    
    % --- C. Calculate Downstream Wall Point (Cancellation) ---
    % Ray from P(n, j) intersects the next wall plate line
    last_x = Px(n_waves, j);
    last_y = Py(n_waves, j);
    
    % Previous Wall Point
    prev_wx = xw(n_waves + j);
    prev_wy = yw(n_waves + j);
    
    % Target Wall Angle (Decreasing: 6, 4, 2, 0)
    wall_theta = (n_waves - j) * theta_step;
    
    % Ray Slope (Right running from last mesh point)
    [~, mu_exit] = GetFlow((n_waves+j)*theta_step, gamma);
    m_ray = tand(((n_waves-j)*theta_step) + mu_exit);
    m_wall = tand(wall_theta);
    
    % Intersect
    wx = (m_ray*last_x - m_wall*prev_wx + prev_wy - last_y) / (m_ray - m_wall);
    wy = m_wall*(wx - prev_wx) + prev_wy;
    
    xw(n_waves + j + 1) = wx;
    yw(n_waves + j + 1) = wy;
end

% --- 4. PLOTTING ---
figure('Color','w', 'Position', [100 100 1200 600]);
hold on; axis equal; grid on;

% Plot Geometry
plot(xw, yw, 'k.-', 'LineWidth', 2); % Upper Wall
yline(0, 'k-.'); % Centerline

% Plot Mesh Lines
% Left-Running (Down)
for i = 1:n_waves
    pX = [xw(i+1)]; pY = [yw(i+1)]; % Start at Wall
    for j = 1:i
        pX(end+1) = Px(i,j); pY(end+1) = Py(i,j);
    end
    pX(end+1) = Cx(i); pY(end+1) = 0; % End at Centerline
    plot(pX, pY, 'Color', [0.6 0.6 0.6]);
end

% Right-Running (Up)
for j = 1:n_waves
    pX = [Cx(j)]; pY = [0]; % Start at Centerline
    for i = j:n_waves
        pX(end+1) = Px(i,j); pY(end+1) = Py(i,j);
    end
    % End at Wall
    pX(end+1) = xw(n_waves + j + 1);
    pY(end+1) = yw(n_waves + j + 1);
    plot(pX, pY, 'Color', [0.6 0.6 0.6]);
end

% --- 5. LABELS (Matching Figure 5 Style) ---
% Throat
text(xw(1)+5, yw(1)/2, '0', 'FontWeight','bold');

% Wall Regions (Simple Expansion)
% Nu = 2, 4, 6, 8
for i = 1:n_waves
    % Centroid of triangle near wall
    if i==1, next_p = Cx(1); else, next_p = Px(i,1); end
    % Roughly between Wall and first mesh point
    mx = (xw(i+1) + next_p)/2 - 2;
    my = (yw(i+1) + Py(i,1))/2 + 2; 
    if i==1, my = (yw(2)+0)/2 + 5; end % adjust first
    
    text(mx, my, num2str(i*2), 'Color', 'b', 'FontSize', 9);
end

% Internal Diamond Regions
% Nu = (i+j)*2
for j = 1:n_waves-1
    for i = j+1:n_waves
        val = (i+j)*2;
        text(Px(i,j), Py(i,j), num2str(val), ...
            'HorizontalAlignment','center', 'FontSize', 8, 'BackgroundColor','w');
    end
end

% Centerline Regions
% Nu = 4, 8, 12, 16
for j = 1:n_waves
    val = (2*j)*2; 
    if j < n_waves
        mx = (Cx(j) + Cx(j+1))/2;
    else
        mx = Cx(end) + 10;
    end
    text(mx, 5, num2str(val), 'Color', 'r', 'FontSize', 9, 'HorizontalAlignment','center');
end

% Cancellation Regions (Near Downstream Wall)
% Nu = 10, 12, 14, 16
for j = 1:n_waves
    val = (n_waves + j)*2; 
    % Rough centroid
    mx = (Px(n_waves,j) + xw(n_waves+j+1))/2;
    my = (Py(n_waves,j) + yw(n_waves+j+1))/2;
    text(mx, my, num2str(val), 'FontSize',8);
end

% Formatting
title('Prandtl-Meyer Function (\nu) Distribution (Kernel Method)');
xlabel('Axial Distance (mm)');
ylabel('Height (mm)');
ylim([-5 70]);

% --- HELPER ---
function [M, mu] = GetFlow(nu_deg, g)
    if nu_deg <= 0.001, M=1; mu=90; return; end
    nu_rad = deg2rad(nu_deg);
    M = fzero(@(m) (sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m^2-1))) - atan(sqrt(m^2-1))) - nu_rad, [1.001 5]);
    mu = rad2deg(asin(1/M));
end