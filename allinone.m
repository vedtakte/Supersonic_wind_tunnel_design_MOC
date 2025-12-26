clear; clc; close all;

% --- 1. DESIGN TARGETS ---
gamma = 1.4;
h_throat = 40.0; % Fixed Throat Half-Height (mm)
M_target = 1.639; % Target Mach Number

% Calculate Target Exit Height based on Area Ratio
AR_target = (1/M_target) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M_target^2))^((gamma+1)/(2*(gamma-1)));
h_exit_target = h_throat * AR_target;

fprintf('DESIGN GOALS:\n');
fprintf('Target Mach: %.4f\n', M_target);
fprintf('Target Area Ratio: %.4f\n', AR_target);
fprintf('Target Exit Height: %.4f mm\n', h_exit_target);
fprintf('------------------------------------------------\n');

% --- 2. SOLVER: FIND REQUIRED L1 ---
% We need to find L1 such that the resulting exit height == h_exit_target
% We define an anonymous function that runs the MOC logic
noc_func = @(L) RunMOC(L, h_throat, gamma) - h_exit_target;

% Find L1 (Guessing between 10mm and 30mm)
options = optimset('TolX', 1e-6);
L1_optimal = fzero(noc_func, [5, 30], options);

fprintf('OPTIMIZED RESULT:\n');
fprintf('Required Length of Plate 1 (AB): %.4f mm\n', L1_optimal);

% --- 3. GENERATE FINAL GEOMETRY ---
[h_final, xw, yw, Cx, Px, Py] = RunMOC(L1_optimal, h_throat, gamma);

fprintf('Final Exit Height: %.4f mm\n', yw(end));
fprintf('Final Area Ratio: %.4f\n', yw(end)/h_throat);
fprintf('------------------------------------------------\n');

% --- 4. GENERATE TABLE DATA ---
theta_walls = [2, 4, 6, 8, 6, 4, 2, 0];
n_segs = 8;
Lengths = zeros(n_segs, 1);
for k=1:n_segs
    Lengths(k) = sqrt((xw(k+1)-xw(k))^2 + (yw(k+1)-yw(k))^2);
end

ResultTable = table((1:8)', Lengths, theta_walls', xw(1:8)', yw(1:8)', xw(2:9)', yw(2:9)', ...
    'VariableNames',{'Seg','Length','Theta','X1','Y1','X2','Y2'});
disp(ResultTable);

% --- 5. PLOT (FIGURE 5 STYLE) ---
figure('Color','w', 'Position', [100 100 1200 600]);
hold on; axis equal; grid on;

% Plot Walls and Centerline
plot(xw, yw, 'k.-', 'LineWidth', 2);
plot(xw, -yw, 'k-', 'LineWidth', 1); % Bottom visual
yline(0, 'k-.');

n_waves = 4;
theta_step = 2;

% Draw Left-Running Waves (Wall -> Centerline)
for i = 1:n_waves
    pX = [xw(i+1)]; pY = [yw(i+1)];
    for j = 1:i
        pX(end+1) = Px(i,j); pY(end+1) = Py(i,j);
    end
    pX(end+1) = Cx(i); pY(end+1) = 0;
    plot(pX, pY, 'Color', [0.6 0.6 0.6]);
end

% Draw Right-Running Waves (Centerline -> Wall)
for j = 1:n_waves
    pX = [Cx(j)]; pY = [0];
    for i = j:n_waves
        pX(end+1) = Px(i,j); pY(end+1) = Py(i,j);
    end
    pX(end+1) = xw(n_waves+j+1); pY(end+1) = yw(n_waves+j+1);
    plot(pX, pY, 'Color', [0.6 0.6 0.6]);
end

% Labels
text(5, h_throat/2, '0', 'FontWeight','bold');
for i = 1:n_waves
    % Wall Regions
    text((xw(i)+xw(i+1))/2, (yw(i)+yw(i+1))/2 - 3, num2str(i*2), 'Color','b', 'FontSize',8);
end
for j = 1:n_waves
    % Centerline Regions (Nu = 4, 8, 12, 16)
    if j<n_waves, mx=(Cx(j)+Cx(j+1))/2; else, mx=Cx(j)+10; end
    text(mx, 2, num2str(j*4), 'Color','r', 'FontSize',8, 'HorizontalAlignment','center');
end
for j=1:n_waves-1
    for i=j+1:n_waves
        text(Px(i,j), Py(i,j), num2str((i+j)*2), 'FontSize',7, 'BackgroundColor','w', 'HorizontalAlignment','center');
    end
end

title(sprintf('Optimized MOC Design (Throat=40mm, M=%.2f)', M_target));
xlabel('Axial Distance (mm)'); ylabel('Height (mm)');
ylim([-5 60]); xlim([0 xw(end)+5]);


% =========================================================================
% CORE MOC SOLVER FUNCTION
% =========================================================================
function [h_exit, xw, yw, Cx, Px, Py] = RunMOC(L1, h_throat, gamma)
    theta_walls = [2, 4, 6, 8, 6, 4, 2, 0];
    n_waves = 4;
    theta_step = 2;
    
    % Init Arrays
    xw = zeros(1, 9); yw = zeros(1, 9);
    Cx = zeros(1, 4);
    Px = zeros(4, 4); Py = zeros(4, 4);
    
    % Point A
    xw(1) = 0; yw(1) = h_throat;
    
    % Point B (Defined by L1)
    xw(2) = xw(1) + L1 * cosd(theta_walls(1));
    yw(2) = yw(1) + L1 * sind(theta_walls(1));
    
    % Calculate remaining geometry using Kernel Logic
    % Note: In Kernel method with variable lengths, we calculate waves first,
    % then intersection determines the next wall point.
    
    % We assume the wall follows the expansion angles initially, 
    % but we need to find WHERE the next wave starts.
    % Actually, for the initial expansion fan (A->E), the wall angles change
    % but the lengths are determined by the reflected characteristics coming back?
    % NO. In Minimum Length Nozzle design:
    % 1. The Expansion contour (A->E) is arbitrary (usually fixed arcs or lengths).
    % 2. The Cancellation contour (E->I) is calculated.
    % To solve for specific Area Ratio with fixed throat, we scale the Expansion portion.
    % So we will scale L2, L3, L4 relative to L1 or keep them equal to L1.
    % Let's assume L1=L2=L3=L4 for the expansion side to keep it smooth.
    
    L_exp = L1; 
    
    % Build Expansion Wall (A->B->C->D->E)
    for i = 1:4
        xw(i+1) = xw(i) + L_exp * cosd(theta_walls(i));
        yw(i+1) = yw(i) + L_exp * sind(theta_walls(i));
    end
    
    % Now Solve Cancellation (Find Cx, Px, Downstream Wall)
    
    % Loop 1: Wall -> Centerline
    for i = 1:n_waves
        wall_idx = i + 1; % Starts at B(2), C(3)...
        sx = xw(wall_idx); sy = yw(wall_idx);
        
        nu_loc = (2*i + 2*(i-1))/2; % Approx flow avg
        if i==1, nu_loc=2; end
        [~, mu] = Flow(nu_loc, gamma);
        lam = (i*2) - mu; % theta - mu
        
        if i==1
            Cx(1) = sx - sy/tand(lam);
        end 
    end
    
    % Loop 2: Internal Grid & Wall Construction
    for j = 1:n_waves
        
        % A. Find C(j)
        if j > 1
            sx = Px(j, j-1); sy = Py(j, j-1);
            % Flow approaching centerline
            nu_loc = (2*j-1)*2; 
            [~, mu] = Flow(nu_loc, gamma);
            lam = 2 - mu; % theta approx 2 deg approaching
            Cx(j) = sx - sy/tand(lam);
        end
        
        % B. Internal Points
        for i = j:n_waves
            if j==1, xu=xw(i+1); yu=yw(i+1); else, xu=Px(i,j-1); yu=Py(i,j-1); end
            if i==j, xl=Cx(j); yl=0; else, xl=Px(i-1,j); yl=Py(i-1,j); end
            
            nu_c = (i+j)*2; theta_c = (i-j)*2;
            [~, mu] = Flow(nu_c, gamma);
            mm = tand(theta_c - mu); mp = tand(theta_c + mu);
            
            xi = (mm*xu - mp*xl + yl - yu)/(mm - mp);
            yi = mm*(xi - xu) + yu;
            Px(i,j) = xi; Py(i,j) = yi;
        end
        
        % C. Find Next Wall Point (Cancellation)
        % Ray from P(4, j) hits wall line
        lx = Px(4, j); ly = Py(4, j);
        
        % Previous Wall Point
        pwx = xw(4+j); pwy = yw(4+j);
        
        % Wall Line Angle
        w_theta = theta_walls(4+j);
        
        % Ray Angle (Right Running)
        nu_ray = (4+j)*2; theta_ray = (4-j)*2;
        [~, mu] = Flow(nu_ray, gamma);
        m_ray = tand(theta_ray + mu);
        m_wall = tand(w_theta);
        
        wx = (m_ray*lx - m_wall*pwx + pwy - ly)/(m_ray - m_wall);
        wy = m_wall*(wx - pwx) + pwy;
        
        xw(4+j+1) = wx;
        yw(4+j+1) = wy;
    end
    
    h_exit = yw(end);
end

function [M, mu] = Flow(nu_deg, g)
    if nu_deg <= 0.001, M=1; mu=90; return; end
    nu_rad = deg2rad(nu_deg);
    try
        M = fzero(@(m) (sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m^2-1))) - atan(sqrt(m^2-1))) - nu_rad, [1.001 5]);
        mu = rad2deg(asin(1/M));
    catch
        M = 1; mu=90;
    end
end