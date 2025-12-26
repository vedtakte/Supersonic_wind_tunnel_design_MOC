clear; clc; close all;

% --- 1. PARAMETERS ---
gamma = 1.4;
L_seg = 20; % Plate length 20mm
theta_step = 2; % 2 degree increments
n_waves = 4; % 4 Expansion waves

% --- 2. GENERATE WALL GEOMETRY (Up to Expansion) ---
h_throat = 39.54; 

% Wall Coords Arrays (Indices: 1=A, 2=B, 3=C...)
n_pts_total = 2 * n_waves + 2; 
xw = zeros(1, n_pts_total);
yw = zeros(1, n_pts_total);

% Point A (Throat)
xw(1) = 0; yw(1) = h_throat;

% Generate Expansion Wall Points (A -> B -> C -> D -> E)
% Angles: 2, 4, 6, 8
current_theta = 0;
for k = 1:n_waves
    current_theta = current_theta + theta_step;
    xw(k+1) = xw(k) + L_seg * cosd(current_theta);
    yw(k+1) = yw(k) + L_seg * sind(current_theta);
end

% --- 3. MOC GRID SOLVER ---
% Px(row, col) -> row = expansion wave index (i), col = reflection index (j)
Px = zeros(n_waves, n_waves);
Py = zeros(n_waves, n_waves);
Cx = zeros(1, n_waves); % Centerline points (y=0)

% We solve Column by Column (Reflected Wave j)
for j = 1:n_waves
    
    % --- STEP A: Find Centerline Point C(j) ---
    % This is where Exp Wave j hits y=0
    
    % Find the point "above" C(j)
    if j == 1
        % For the first wave, the source is Wall Point B (Index 2)
        src_x = xw(2); 
        src_y = yw(2);
        
        % Flow in this region is simple expansion nu=2, theta=2
        theta_loc = 1 * theta_step;
        nu_loc = theta_loc;
    else
        % For subsequent waves, source is P(j, j-1)
        src_x = Px(j, j-1);
        src_y = Py(j, j-1);
        
        % Flow approaching centerline
        % Region is downstream of P(j, j-1)
        % Theta approx theta_step, Nu approx (2*j)*theta_step
        theta_loc = theta_step;
        nu_loc = (2*j - 1) * theta_step; 
    end
    
    [~, mu] = GetFlow(nu_loc, gamma);
    
    % Characteristic slope (Left-Running: theta - mu)
    lambda = theta_loc - mu;
    
    % Solve intersection with y=0
    % 0 - src_y = m * (x - src_x)  ->  x = src_x - src_y/m
    Cx(j) = src_x - src_y / tand(lambda);
    
    
    % --- STEP B: March Upwards (Reflected Wave j) ---
    for i = j:n_waves
        % Find Intersection P(i,j)
        
        % 1. Get Lower-Left Source (Right-Running Stream)
        if i == j
            xl = Cx(j); yl = 0; % Comes from Centerline
        else
            xl = Px(i-1, j); yl = Py(i-1, j);
        end
        
        % 2. Get Upper-Left Source (Left-Running Stream)
        if j == 1
            xu = xw(i+1); yu = yw(i+1); % Comes from Wall
        else
            xu = Px(i, j-1); yu = Py(i, j-1);
        end
        
        % 3. Properties in the target cell (i,j)
        theta_cell = (i - j) * theta_step;
        nu_cell = (i + j) * theta_step;
        [~, mu_cell] = GetFlow(nu_cell, gamma);
        
        % 4. Intersect Lines
        m_minus = tand(theta_cell - mu_cell); % Slope from Upper
        m_plus  = tand(theta_cell + mu_cell); % Slope from Lower
        
        x_int = (m_minus*xu - m_plus*xl + yl - yu) / (m_minus - m_plus);
        y_int = m_minus*(x_int - xu) + yu;
        
        Px(i, j) = x_int;
        Py(i, j) = y_int;
    end
    
    % --- STEP C: Wall Cancellation (Find next wall point) ---
    % Connect P(n, j) to Wall
    last_p_x = Px(n_waves, j);
    last_p_y = Py(n_waves, j);
    
    % Target Wall Angle
    wall_theta = (n_waves - j) * theta_step; 
    
    % Previous Wall Point
    prev_wx = xw(n_waves + j);
    prev_wy = yw(n_waves + j);
    
    % Right running wave slope
    m_wave = tand(((n_waves - j)*theta_step) + mu_cell); 
    m_wall = tand(wall_theta);
    
    % Intersection
    wx = (m_wave*last_p_x - m_wall*prev_wx + prev_wy - last_p_y) / (m_wave - m_wall);
    wy = m_wall*(wx - prev_wx) + prev_wy;
    
    xw(n_waves + j + 1) = wx;
    yw(n_waves + j + 1) = wy;
end

% --- 4. PLOTTING ---
figure('Color','w', 'Position', [100 100 1200 500]);
hold on; axis equal; grid on;

% Plot Walls
plot(xw(1:end-1), yw(1:end-1), 'k.-', 'LineWidth', 2);
yline(0, 'k-.');

% --- PLOT LEFT-RUNNING WAVES (Expansion) ---
% These go: Wall -> Mesh -> Centerline
for i = 1:n_waves
    % Start at Wall
    pathX = xw(i+1); 
    pathY = yw(i+1);
    
    % Go through mesh points
    for j = 1:i
        pathX(end+1) = Px(i, j);
        pathY(end+1) = Py(i, j);
    end
    
    % *** FIX: CONNECT TO CENTERLINE ***
    pathX(end+1) = Cx(i);
    pathY(end+1) = 0;
    
    plot(pathX, pathY, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

% --- PLOT RIGHT-RUNNING WAVES (Reflection) ---
% These go: Centerline -> Mesh -> Wall
for j = 1:n_waves
    % Start at Centerline
    pathX = Cx(j);
    pathY = 0;
    
    % Go through mesh points
    for i = j:n_waves
        pathX(end+1) = Px(i, j);
        pathY(end+1) = Py(i, j);
    end
    
    % Connect to Wall
    pathX(end+1) = xw(n_waves + j + 1);
    pathY(end+1) = yw(n_waves + j + 1);
    
    plot(pathX, pathY, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

% --- 5. LABELS ---
text(5, 20, '0', 'FontWeight','bold'); % Region 0

% Wall Region Labels
for i=1:n_waves
    % Approximate position: Midpoint of wall segment
    xm = (xw(i)+xw(i+1))/2; ym = (yw(i)+yw(i+1))/2;
    % nu = 2, 4, 6, 8
    % text(xm, ym-5, num2str(i*2), 'Color', 'b', 'FontSize', 8); 
end

% Centerline Labels
for j = 1:n_waves
    % nu = 2, 6, 10, 14 (As per your PDF Figure 5)
    % or 0, 4, 8, 12 depending on specific convention.
    % Based on your PDF, centerline values are 4, 8, 12, 16
    val = 2 * j * theta_step;
    
    % Position text between C(j) and C(j+1)
    if j < n_waves
        xm = (Cx(j) + Cx(j+1))/2;
    else
        xm = Cx(end) + 15;
    end
    text(xm, 2, num2str(val), 'Color', 'r', 'FontSize', 9, 'HorizontalAlignment', 'center');
end

% Internal Cell Labels
for j = 1:n_waves-1
    for i = j+1:n_waves
        val = (i+j)*theta_step;
        text(Px(i,j), Py(i,j), num2str(val), 'HorizontalAlignment', 'center', 'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);
    end
end

title('Figure 5: Prandtl-Meyer Function (\nu) in each cell');
xlabel('Axial Distance (mm)');
ylabel('Height (mm)');
xlim([0 165]); ylim([0 55]);

% --- HELPER FUNCTIONS ---
function [M, mu] = GetFlow(nu_deg, g)
    if nu_deg <= 0.001, M=1; mu=90; return; end
    nu_rad = deg2rad(nu_deg);
    M = fzero(@(m) (sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m^2-1))) - atan(sqrt(m^2-1))) - nu_rad, [1.001 5]);
    mu = rad2deg(asin(1/M));
end