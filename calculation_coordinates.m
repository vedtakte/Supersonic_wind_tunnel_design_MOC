clear; clc;

% --- 1. INPUTS ---
gamma = 1.4;
L_plate = 20; % mm (Plate Length)
theta_schedule = [2, 4, 6, 8, 6, 4, 2, 0]; % Wall angles in degrees

% --- 2. AERODYNAMIC SIZING ---
% Calculate Max Turn (half of total expansion)
theta_max = max(theta_schedule);
nu_exit = 2 * deg2rad(theta_max);

% Solve for Mach Number at Exit
M_exit = fzero(@(M) PM_Function(M, gamma) - nu_exit, [1.01, 3.0]);

% Calculate Area Ratio (A_exit / A_throat)
AR = (1/M_exit) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M_exit^2))^((gamma+1)/(2*(gamma-1)));

% --- 3. GEOMETRIC INTEGRATION ---
% We need to find the specific Throat Height (y0) that enables this Area Ratio
% given that our Wall Lengths and Angles are fixed.
% y_exit = y_throat + sum(L * sin(theta))
% AR = y_exit / y_throat
% -> y_throat = sum(L * sin(theta)) / (AR - 1)

delta_y_total = sum(L_plate * sind(theta_schedule));
h_throat = delta_y_total / (AR - 1);
h_exit = h_throat * AR;

% --- 4. COORDINATE CALCULATION LOOP ---
num_plates = length(theta_schedule);
Plate_Names = cell(num_plates, 1);
Lengths = zeros(num_plates, 1);
Angles = zeros(num_plates, 1);
Start_X = zeros(num_plates, 1);
Start_Y = zeros(num_plates, 1);
End_X = zeros(num_plates, 1);
End_Y = zeros(num_plates, 1);

% Initialize Start Point (Point A)
current_x = 0;
current_y = h_throat;

% Generate Alphabet for Point Names (A, B, C...)
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

for i = 1:num_plates
    % Define Segment Name (e.g., AB, BC)
    p1 = alphabet(i);
    p2 = alphabet(i+1);
    Plate_Names{i} = [p1 p2];
    
    % Store Properties
    Lengths(i) = L_plate;
    Angles(i) = theta_schedule(i);
    
    % Store Start Coordinates
    Start_X(i) = current_x;
    Start_Y(i) = current_y;
    
    % Calculate End Coordinates
    current_x = current_x + L_plate * cosd(theta_schedule(i));
    current_y = current_y + L_plate * sind(theta_schedule(i));
    
    % Store End Coordinates
    End_X(i) = current_x;
    End_Y(i) = current_y;
end

% --- 5. CREATE TABLE ---
GeomTable = table(Plate_Names, Lengths, Angles, Start_X, Start_Y, End_X, End_Y, ...
    'VariableNames', {'Plate', 'Length_mm', 'Theta_deg', 'X1', 'Y1', 'X2', 'Y2'});

% Display Main Table
fprintf('\nTable 3: Lengths of Plates, Angles, and Coordinates\n');
disp(GeomTable);

% Display Critical Dimensions (Like OA and IP in your PDF)
fprintf('--------------------------------------------------\n');
fprintf('CRITICAL DIMENSIONS (Vertical Heights)\n');
fprintf('--------------------------------------------------\n');
fprintf('Throat Half-Height (OA):  %.3f mm\n', h_throat);
fprintf('Exit Half-Height   (IP):  %.3f mm\n', h_exit);
fprintf('Total Axial Length:       %.3f mm\n', current_x);
fprintf('--------------------------------------------------\n');

% --- Helper Function ---
function nu = PM_Function(M, g)
    term1 = sqrt((g+1)/(g-1));
    term2 = atan(sqrt((g-1)/(g+1) * (M^2 - 1)));
    term3 = atan(sqrt(M^2 - 1));
    nu = term1 * term2 - term3;
end