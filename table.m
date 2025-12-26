clear; clc;

% --- 1. Setup Parameters ---
gamma = 1.4;
step_size = 2;
nu_max = 16;

% Initialize Arrays for Calculation
nu_vals = 0:step_size:nu_max;
mu_vals = zeros(size(nu_vals));
M_vals = zeros(size(nu_vals));

% --- 2. Calculate State Properties (nu, mu, M) ---
for i = 1:length(nu_vals)
    nu_deg = nu_vals(i);
    nu_rad = deg2rad(nu_deg);
    
    % Solve for Mach (M)
    if nu_deg == 0
        M = 1.0;
    else
        % Invert P-M function to find Mach
        M = fzero(@(m) PM_Function(m, gamma) - nu_rad, [1.0001, 5.0]);
    end
    M_vals(i) = M;
    
    % Calculate Mach Angle (mu) in degrees
    mu_vals(i) = rad2deg(asin(1/M));
end

% --- 3. Calculate Characteristic Angles and Store in Vectors ---
% We pre-allocate vectors to store the results for the table
num_rows = length(nu_vals) - 1;
Nu1_vec = zeros(num_rows, 1);
Mu1_vec = zeros(num_rows, 1);
Nu2_vec = zeros(num_rows, 1);
Mu2_vec = zeros(num_rows, 1);
Lambda_vec = zeros(num_rows, 1);

for i = 1:length(nu_vals)-1
    % State 1 (Previous)
    Nu1_vec(i) = nu_vals(i);
    Mu1_vec(i) = mu_vals(i);
    
    % State 2 (Next)
    Nu2_vec(i) = nu_vals(i+1);
    Mu2_vec(i) = mu_vals(i+1);
    
    % Apply the averaging formula: Lambda = 0.5 * [ (mu1 - nu1) + (mu2 - nu2) ]
    term1 = Mu1_vec(i) - Nu1_vec(i);
    term2 = Mu2_vec(i) - Nu2_vec(i);
    Lambda_vec(i) = 0.5 * (term1 + term2);
end

% --- 4. Create and Display the Table ---
CharacteristicTable = table(Nu1_vec, Mu1_vec, Nu2_vec, Mu2_vec, Lambda_vec, ...
    'VariableNames', {'nu1_deg', 'mu1_deg', 'nu2_deg', 'mu2_deg', 'Lambda_char_deg'});

% Display the table in Command Window
disp('Characteristic Angles Table:');
disp(CharacteristicTable);

% --- Helper Function ---
function nu = PM_Function(M, g)
    term1 = sqrt((g+1)/(g-1));
    term2 = atan(sqrt((g-1)/(g+1) * (M^2 - 1)));
    term3 = atan(sqrt(M^2 - 1));
    nu = term1 * term2 - term3;
end