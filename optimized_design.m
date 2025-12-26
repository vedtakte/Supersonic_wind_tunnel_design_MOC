clear; clc; close all;

% --- 1. DESIGN PARAMETERS (OPTIMIZED) ---
gamma = 1.4;
h_throat = 40.00;       % Throat Half-Height (mm)
L_exp = 14.35;          % Optimized Plate Length for AB, BC, CD, DE
theta_step = 2;         % Angle step (deg)
n_waves = 4;            % Number of expansion waves

% --- 2. GENERATE EXPANSION WALL (A -> E) ---
xw = zeros(1, 9); yw = zeros(1, 9);
xw(1) = 0; yw(1) = h_throat; % Point A

% Build points B, C, D, E
current_theta = 0;
for i = 1:n_waves
    current_theta = current_theta + theta_step;
    xw(i+1) = xw(i) + L_exp * cosd(current_theta);
    yw(i+1) = yw(i) + L_exp * sind(current_theta);
end

% --- 3. MOC SOLVER (Method of Characteristics) ---
Px = zeros(n_waves, n_waves); Py = zeros(n_waves, n_waves);
Cx = zeros(1, n_waves); % Centerline

for j = 1:n_waves
    % A. Find Centerline Intersections (C)
    if j == 1
        sx = xw(2); sy = yw(2); nu_loc = 2; theta_loc = 2;
    else
        sx = Px(j, j-1); sy = Py(j, j-1);
        nu_loc = (2*j-1)*2; theta_loc = 2;
    end
    [~, mu] = GetFlow(nu_loc, gamma);
    Cx(j) = sx - sy / tand(theta_loc - mu);

    % B. Internal Grid Points (P)
    for i = j:n_waves
        if j==1, xu=xw(i+1); yu=yw(i+1); else, xu=Px(i,j-1); yu=Py(i,j-1); end
        if i==j, xl=Cx(j); yl=0; else, xl=Px(i-1,j); yl=Py(i-1,j); end
        
        nu_c = (i+j)*2; theta_c = (i-j)*2;
        [~, mu] = GetFlow(nu_c, gamma);
        mm = tand(theta_c - mu); mp = tand(theta_c + mu);
        
        Px(i,j) = (mm*xu - mp*xl + yl - yu)/(mm - mp);
        Py(i,j) = mm*(Px(i,j) - xu) + yu;
    end

    % C. Cancellation Wall Points (F, G, H, I)
    lx = Px(n_waves, j); ly = Py(n_waves, j); % Last mesh point
    pwx = xw(n_waves+j); pwy = yw(n_waves+j); % Previous wall point
    
    w_theta = (n_waves - j) * theta_step; % Target wall angle
    [~, mu_r] = GetFlow((n_waves+j)*2, gamma);
    theta_r = (n_waves-j)*2;
    
    m_ray = tand(theta_r + mu_r);
    m_wall = tand(w_theta);
    
    xw(n_waves+j+1) = (m_ray*lx - m_wall*pwx + pwy - ly)/(m_ray - m_wall);
    yw(n_waves+j+1) = m_wall*(xw(n_waves+j+1) - pwx) + pwy;
end

% Add Exit Extension
xw = [xw, xw(end)+30]; yw = [yw, yw(end)];

% --- 4. PLOTTING ---
figure('Color','w', 'Position', [100 100 1000 500]);
hold on; axis equal; grid on;

% Plot Walls
plot(xw, yw, 'k.-', 'LineWidth', 2, 'MarkerSize', 12);       % Upper
plot(xw, -yw, 'k.-', 'LineWidth', 2, 'MarkerSize', 12);      % Lower
yline(0, 'k-.', 'LineWidth', 1);                             % Centerline

% Labels
text(-5, h_throat, 'Throat', 'HorizontalAlignment','right');
text(xw(end), yw(end)+2, 'Test Section', 'HorizontalAlignment','center');
labels = {'A','B','C','D','E','F','G','H','I'};
for k=1:9
    text(xw(k), yw(k)+3, labels{k}, 'FontWeight','bold', 'HorizontalAlignment','center');
end

% Annotations
dim_str = sprintf('Exit Height: %.2f mm\nArea Ratio: %.3f', yw(end)*2, yw(end)/h_throat);
text(xw(end)-10, 10, dim_str, 'BackgroundColor','w', 'EdgeColor','k');

title('Supersonic Wind Tunnel Nozzle Design (M=1.64, Throat=40mm)');
xlabel('Axial Length (mm)'); ylabel('Height (mm)');
ylim([-65 65]); xlim([-10 max(xw)+10]);

% --- HELPER ---
function [M, mu] = GetFlow(nu_deg, g)
    if nu_deg <= 0.001, M=1; mu=90; return; end
    nu_rad = deg2rad(nu_deg);
    M = fzero(@(m) (sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m^2-1))) - atan(sqrt(m^2-1))) - nu_rad, [1.001 5]);
    mu = rad2deg(asin(1/M));
end