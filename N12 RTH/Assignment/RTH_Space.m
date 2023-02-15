% This Octave/MATLAB Script Performs Thermohydraulic Calculations
% For a NERVA Type NTP in Space

clc;
clear;
close all;

% Define Reactor Core Data
D_ele = 1.905E-2; % m
s_ele = 10.8E-3; % m
D_hole = 2E-3; % m
D_sep = 25E-3 / 6; % m
t_clad = 0.25E-3; % m
D_chnl = D_hole - 2 * t_clad; % m
A_ele = 3 * s_ele * D_ele / 2; % m^2

Q_c = 2.25E9; % W
Q_ele = 0.65E6; % W
H_ele = 0.89; % m
H_ex = 0.93463; % m
Num_ele = ceil(Q_c / Q_ele * 7/6); % Includes Tie Tubes in a 6 : 1 Ratio
m_UC = Num_ele * 13.63E3 * (A_ele - 19 * 0.25 * pi * D_chnl.^2) * H_ele; % Kg

T_mar = 40; % K
SF_PV = 1.5;
Sheath_W = 1; % Units of t_PV

P_init = 4.7338E6; % Pa
T_init = 2642.0; % K
m_dot = 32; % Kg / s

% Define De Laval Nozzle Geometry
C_deg = deg2rad(35); % rad
D_deg = deg2rad(15); % rad
Ae_by_A_t = 100;

% Define Turbine Data
Tmax_turbin = 1373; % K (Maximum Uncooled Turbine Blade Limit)
PRATIO = 1.61^3; % Raised to the Number of Stages
M_y = 6.5E-3; % Units of m_dot

% Define H2 Material Data Functions
% Load H2 Material Data
Cp_H2_dat = csvread("H2_Cp_5MPa.csv");
mu_H2_dat = csvread("H2_mu_5MPa.csv");
k_H2_dat = csvread("H2_k_5MPa.csv");
Pr_H2_dat = csvread("H2_Pr_5MPa.csv");

% Isobaric Specific Heat vs Temperature
Cp_H2 = @ (T) interp1(Cp_H2_dat(:, 1), Cp_H2_dat(:, 2), T, "cubic"); % J / Kg * K

% Dynamic Viscosity vs Temperature
mu_H2 = @ (T) interp1(mu_H2_dat(:, 1), mu_H2_dat(:, 2), T, "cubic"); % Pa * s

% Thermal Conductivity vs Temperature
k_H2 = @ (T) interp1(k_H2_dat(:, 1), k_H2_dat(:, 2), T, "cubic"); % W / m * K

% Prandtl Number vs Temperature
Pr_H2 = @ (T) interp1(Pr_H2_dat(:, 1), Pr_H2_dat(:, 2), T, "cubic");

% Ratio of Specific Heats
gamma_H2 = 7 / 5; % Ideal Diatomic Gas (Actual Varies 1.597 to 1.318)

% Specific Gas Constant
R_H2 = 8.314462 / 2E-3; % J / Kg * K


% Define MgO Material Data Functions
% Load MgO Material Data
k_MgO_dat = csvread("MgO_k.csv");

% Thermal Conductivity vs Temperature
k_MgO = @ (T) interp1(k_MgO_dat(:, 1), k_MgO_dat(:, 2), T, "cubic", "extrap"); % W / m * K

% Maximum Working Temperature
Tmax_MgO = 3125 - T_mar; % K


% Define Al2O3 Material Data Functions
% Load Al2O3 Material Data
k_Al2O3_dat = csvread("Al2O3_k.csv");

% Thermal Conductivity vs Temperature
k_Al2O3_dat = @ (T) interp1(k_Al2O3_dat(:, 1), k_Al2O3_dat(:, 2), T, "cubic", "extrap"); % W / m * K

% Maximum Working Temperature
Tmax_Al2O3 = 2273 - T_mar; % K


% Define UC Material Data Functions
% Load UC Material Data
k_UC_dat = csvread("UC_k.csv");

% Thermal Conductivity vs Temperature
k_UC = @ (T) interp1(k_UC_dat(:, 1), k_UC_dat(:, 2), T, "cubic", "extrap"); % W / m * K

% Maximum Working Temperature
Tmax_UC = 2780 - T_mar; % K

% Define A-286 Steel Mechanical Data
% Allowable Stress at 645 K for S66286 of ASTM SA-638
SA_A286 = 241E6; % Pa


% Define Helper Function for Circular Cross Sectional Area
A_circ = @ (r1, r2) (pi .* (r2.^2 - r1.^2)); % m^2

% Define Helper Function for Reynold's Number
f_Re = @ (M, R1, R2,  T) (2 .* M ./ (pi .* (R1 + R2) .* mu_H2(T)));

% Define Helper Function for Darcy Factor
f_dw = @ (Re) (((64 ./ Re) .* (Re < 3E3)) + ((0.316 .* Re.^-0.25) .* ((Re > 3E3) & (Re < 3E4))) + ((0.184 .* Re.^-0.2) .* (Re > 3E4)));

% Define Helper Function for Ideal Gas Density
rho = @ (P, R, T) (P ./ (R .* T)); % Kg / m^3

% Define Helper Function for Isentropic Nozzle Area Ratios
A_by_A_t = @ (M, g) ((1 ./ M) .* ((g + 1) ./ (2 + (g - 1) .* M.^2)).^((g + 1) ./ (-2 .* g + 2)));

% Define Helper Function for Isentropic Nozzle Subsonic Mach Number
M_A_sub = @(r, g) (fzero(@ (x) (A_by_A_t (x, g) - r), [0, 1]));

% Define Helper Function for Isentropic Nozzle Supersonic Mach Number
M_A_sup = @(r, g) (fzero(@ (x) (A_by_A_t (x, g) - r), [1, 100]));

% Define Helper Function for Isentropic Total Temperature
T0_by_T = @ (M, g) (1 + (g - 1) .* 0.5 .* M.^2);


% Initialize PMAT Matrix Entries 1-7
% Records Pressure - Mass Flow - Area - Temperature at Various Points in Reactor
% 1 - Chamber
% 2 - LH2 Tank
% 3 - Pump Outlet / Sheath Entry
% 4 - Sheath Exit / Reactor Plenum
% 5 - Turbine Inlet
% 6 - Turbine Outlet
% 7 - Nozzle Outlet
PMAT = zeros(7, 4);

% Setup Point 1 - Chamber
PMAT(1, 1) = P_init; % Pa
PMAT(1, 2) = m_dot; % Kg / s
PMAT(1, 3) = Num_ele * A_ele; % m^2 (6 : 1 Tie Tube Ratio)
PMAT(1, 4) = T_init; % K

% Calculate Nozzle Inlet Conditions
R_in = sqrt(PMAT(1, 3) / pi); % m
M_in = PMAT(1, 2) / (PMAT(1, 3) * rho(PMAT(1, 1), R_H2, PMAT(1, 4)) * sqrt(gamma_H2 * R_H2 * PMAT(1, 4)));
T0 = PMAT(1, 4) * T0_by_T(M_in, gamma_H2); % K
P0 = PMAT(1, 1) * T0_by_T(M_in, gamma_H2)^(gamma_H2 / (gamma_H2 - 1)); % Pa

% Calculate Nozzle Throat Data
A_t = PMAT(1, 3) / A_by_A_t(M_in, gamma_H2); % m^2
R_t = sqrt(A_t / pi); % m

% Calculate Nozzle Exit Data
Ae = A_t * Ae_by_A_t; % m^2
R_e = sqrt(Ae / pi); % m
M_e = M_A_sup(Ae_by_A_t, gamma_H2); % m

% Calculate Nozzle Section Heights
Z_C = (R_in - R_t) / tan(C_deg); % m (Converging Section)
Z_D = (R_e - R_t) / tan(D_deg); % m (Diverging Section)

% Calculate Diverging Side Radiative Separation Point
A_sep = A_t * A_by_A_t(fzero(@ (x) (T0 / T0_by_T(x, gamma_H2) - Tmax_Al2O3), [0 100]), gamma_H2); % m^2
R_sep = sqrt(A_sep / pi); % m
Z_sep = (R_sep - R_t) / tan(C_deg); % m


% Setup Point 2 - LH2 Tank
PMAT(2, 1) = 0.225E6; % Pa
PMAT(2, 2) = PMAT(1, 2); % Kg / s
PMAT(2, 3) = A_circ(0, 102.26E-3); % m^2 (NPS 4 Sch 40 Pipe)
PMAT(2, 4) = 20; % K


% Setup Point 3 - Post LH2 Pump - Sheath Entry
PMAT(3, 1) = 5E6; % Pa
PMAT(3, 2) = PMAT(1, 2); % Kg / s

% Calculate Chamber Vessel Thickness and Cooling Sheath Dimensions
t_PV = (SF_PV * PMAT(3, 1) * R_in) / (SA_A286 + 0.4 * PMAT(3, 1) * SF_PV); % m
t_PV = ceil(t_PV * 1E3) / 1E3; % Round to Nearest 1E-3 m

% Finish Point 3 Setup
PMAT(3, 3) = A_circ(R_t + t_PV, R_t + (1 + Sheath_W) * t_PV); % m^2
PMAT(3, 4) = 20; % K


% Calculate Pump Input Work at 20 K and 5 MPa
PWR_Pmp = (PMAT(1, 2) * 1 / 76.2761 * (PMAT(3, 1) - PMAT(2, 1))) / (0.9 * 0.7258); % W


% Calculate Regeneratively Cooled Nozzle Data Using Conduction Approximation
Q_rgen_cond = 15 * pi * (R_in + R_t) * Z_C / cos(C_deg) * T0 / t_PV; % W (A-286 Conductivity is 15 W / m * K)
T_rgen_cond = fzero(@ (x) (integral(Cp_H2, PMAT(3, 4), x, "AbsTol", 1E-15) - (Q_rgen_cond / PMAT(1, 2))), [PMAT(3, 4) PMAT(1, 4)]); % K


% Calculate Regenratively Cooled Nozzle Using Bartz Heat Transfer Approximation
% Setup Helper Functions
R_rgen = @ (z) (R_t + tan(C_deg) .* z); % m
A_rgen = @ (z) A_circ(0, R_rgen(z)); % m^2
T_rgen_z = @ (z) (T0 ./ T0_by_T(M_A_sub(A_rgen(z) ./ A_t, gamma_H2), gamma_H2)); % K
dh_f_rgen = @ (z) (0.026 .* (k_H2(T_rgen_z(z)) ./ (2 .* R_rgen(z))) .* (f_Re(PMAT(1, 2), 0, R_rgen(z), T_rgen_z(z)).^0.8) .* (Pr_H2(T_rgen_z(z)).^0.4)); % W / m^2 * K

% Integrate Heat Transfer Coefficient Over Nozzle Converging Section Area
Q_rgen_bartz = 0; % Pa
div = 1024;
for i = 1:div
  Q_rgen_bartz += dh_f_rgen(i * (Z_C / div)) * 2 * pi * R_rgen(i * (Z_C / div)) * (tan(C_deg) / cos(C_deg)) * T_rgen_z(i * (Z_C / div)) * (Z_C / div);
endfor

% Calculate Cooling Sheath Outlet Temperature
T_rgen_bartz = fzero(@ (x) (integral(Cp_H2, PMAT(3, 4), x, "AbsTol", 1E-15) - (Q_rgen_bartz / PMAT(1, 2))), [PMAT(3, 4) PMAT(1, 4)]); % K

% Select Maximum of Sheath Outlet Temperature of Two Methods as Outlet
T_rgen = max(T_rgen_cond, T_rgen_bartz); % K

% Setup and Calculate Maximum Pressure Loss in Nozzle Sheath
##R_rgen = @ (z) (R_t + tan(C_deg) .* z); % m
##A_rgen = @ (z) A_circ(R_rgen(z), (R_rgen(z) + (1 + Sheath_W) .* t_PV)); % m^2
##T_rgen_z = @ (z) (PMAT(3, 4) + ((T_rgen_cond - PMAT(3, 4)) ./ Z_C) .* z); % K

##S_w_rgen = @ (z) (2 .* pi .* (2 .* R_rgen(z) + (1 + Sheath_W) .* t_PV)); % m
##dP_rgen = @ (z) ((f_dw(f_Re(PMAT(1, 2), 0, R_rgen(z), T_rgen_z(z))) .* ((PMAT(1, 2) ./ A_rgen(z)).^2)) ./ (8 .* rho(PMAT(3, 1), R_H2, T_rgen_z(z)) .* (A_rgen(z) ./ S_w_rgen(z)))); % Pa
##P_loss = 0; % Pa
##dP_rgen = @ (z) ((PMAT(3, 2) .^ 2) .* (1 ./ (A_rgen(z) .* rho(PMAT(3, 1) - P_loss, R_H2, T_rgen_z(z))))); % Pa
##for i = 0:2048
##  P_loss += dP_rgen(i * (Z_C / (cos(C_deg) * 2048))) + 0 * (Z_C / (cos(C_deg) * 2048));
##endfor
DP_rgen = fzero(@ (x) (x - (((PMAT(3, 2) ./ PMAT(3, 3)).^2) .* (rho(PMAT(3, 1) - x, R_H2, T_rgen).^-1 - rho(PMAT(3, 1), R_H2, PMAT(3, 4)).^-1))), [0 1E6]); % Pa

% Setup Point 4 - Reactor Core Plenum
PMAT(4, 1) = PMAT(3, 1) - DP_rgen; % Pa
PMAT(4, 2) = PMAT(1, 2); % Kg / s
PMAT(4, 3) = A_circ(R_in + t_PV, R_in + (1 + Sheath_W) * t_PV); % m^2
PMAT(4, 4) = T_rgen; % K


% Setup Reactor Core Calculation
% Initialize Reactor Temperature and Power Profile Container
div = 1024;
TBCSM = zeros(div, 4);
DP_c = 0; % Pa

% Calculate Reactor Element Channel Data and Setup Helper Functions
A_chnl = A_circ(0, D_chnl / 2); % m^2
A_flow = 19 * A_chnl; % m^2
M_chnl = (PMAT(1, 2) * Q_ele) / (19 * Q_c); % Kg / s
Q_fc = @ (z) (H_ele ./ div) .* ((Q_ele ./ 19) .* cos(pi .* ((H_ele / 2) - z) ./ H_ex)); % W
HTC  = @ (T) ((k_H2(T) ./ D_chnl) .* 0.205 .* (f_Re(M_chnl, 0, D_chnl / 2, T).^0.8) .* (Pr_H2(T).^0.4)); % W / m^2 * K

% Calculate Inlet Initial Temperatures
TBCSM(1, 1) = PMAT(4, 4); % K
TBCSM(1, 2) = fzero(@ (x) (TBCSM(1, 1) - x + (Q_fc(H_ele ./ div) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(1, 1)) .* ((x ./ TBCSM(1, 1)).^-0.55))))), [PMAT(3, 4) Tmax_MgO]); % K
##TBCSM(1, 2) = (TBCSM(1, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(1, 1))))); % K
TBCSM(1, 3) = ((Q_fc(H_ele ./ div) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(1, 2)))) + TBCSM(1, 2); % K
TBCSM(1, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(1, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(1, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(1, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(1, 3)); % K
TBCSM(2, 1) = TBCSM(1, 1) + Q_fc(H_ele ./ div) ./ (M_chnl .* Cp_H2(TBCSM(1, 1))); % K
DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(1, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(1, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa

% Calculate Full TBCSM Matrix
for i = 2:(div - 1)
  TBCSM(i, 2) = fzero(@ (x) (TBCSM(i, 1) - x + (Q_fc(i .* H_ele ./ div) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(i, 1)) .* ((x ./ TBCSM(i, 1)).^-0.55))))), [PMAT(3, 4) Tmax_MgO]); % K
##  TBCSM(i, 2) = (TBCSM(i, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(i, 1))))); % K
  TBCSM(i, 3) = ((Q_fc(i .* H_ele ./ div) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(i, 2)))) + TBCSM(i, 2); % K
  TBCSM(i, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(i, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(i, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(i, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(i, 3)); % K
  TBCSM(i + 1, 1) = TBCSM(i, 1) + Q_fc(i .* H_ele ./ div) ./ (M_chnl .* Cp_H2(TBCSM(i, 1))); % K
  DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(i, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(i, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa
endfor

% Calculate Core Outlet Values for TBCSM Matrix
TBCSM(div, 2) = fzero(@ (x) (TBCSM(div, 1) - x + (Q_fc(H_ele) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(div, 1)) .* ((x ./ TBCSM(div, 1)).^-0.55))))), [PMAT(3, 4) Tmax_MgO]); % K
##TBCSM(div, 2) = (TBCSM(div, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(div, 1))))); % K
TBCSM(div, 3) = ((Q_fc(H_ele) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(div, 2)))) + TBCSM(div, 2); % K
TBCSM(div, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(div, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(div, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(div, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(div, 3)); % K
DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(div, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(div, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa

% Find Maximum Temperatures Locations
[TmaxM, TmaxInd] = max(TBCSM);
TmaxNorm = TmaxM ./ [Tmax_turbin Tmax_MgO Tmax_UC Tmax_UC];
TmaxZ = (H_ele / 2) .- TmaxInd .* (H_ele / div); % m

% Calculate Core Pressure Loss
##P_loss = 0; % Pa
##dP_rgen = @ (z) ((M_chnl .^ 2) .* (1 ./ (A_chnl .* rho(PMAT(4, 1) - P_loss, R_H2, interp1(0:(H_ele / div):(H_ele - (H_ele / div)), TBCSM(:, 1), z, "cubic"))))); % Pa
##dP_rgen = @ (z, i) (M_chnl .* (M_chnl ./ (A_chnl.^2 .* rho(PMAT(4, 1) - P_loss, R_H2, TBCSM(i, 1))))); % Pa
##rho_c = @ (z) (rho(PMAT(4, 1), R_H2, TBCSM(i, 1)));
##u_c = @ (z) (M_chnl ./ (A_chnl .* rho_c(i)));
##for i = 1:div
##  P_loss += (rho_c(i) .* u_c(i)) .* ((u_c(div) - u_c(1)) / div);
##endfor
DP_c = fzero(@ (x) (x - (((M_chnl ./ A_chnl).^2) .* (rho(PMAT(4, 1) - x, R_H2, TBCSM(1024, 1)).^-1 - rho(PMAT(4, 1), R_H2, PMAT(4, 4)).^-1))), [0 1E6]); % Pa

% Calculate Mass Flow Tap from LH2 Tank Required to Cool Turbine Inlet Gas to Turbine Inlet Temperature
m_turb_cool = (M_y .* PMAT(4, 2) .* integral(Cp_H2, Tmax_turbin, PMAT(1, 4), "AbsTol", 1E-15)) ./ integral(Cp_H2, PMAT(2, 4), Tmax_turbin, "AbsTol", 1E-15); % Kg / s

% Setup Point 5 - Turbine Inlet
PMAT(5, 1) = PMAT(4, 1) - DP_c; % Pa
PMAT(5, 2) = M_y * PMAT(4, 2) + m_turb_cool; % Kg / s
PMAT(5, 3) = A_circ(0, 102.26E-3); % m^2 (NPS 4 Sch 40 Pipe)
PMAT(5, 4) = Tmax_turbin; % K

% Calculate Turbine Output and Net Eletrical Power
PWR_Turb = 0.9 .* 0.7171 .* PMAT(5, 2) .* Cp_H2(PMAT(5, 4)) .* PMAT(5, 4) .* ((PRATIO.^((gamma_H2 - 1) ./ gamma_H2)) - 1); % W
PWR_net = PWR_Turb - PWR_Pmp; % W

% Setup Point 6 - Turbine Outlet
PMAT(6, 1) = PMAT(5, 1) / PRATIO; % Pa
PMAT(6, 2) = M_y * PMAT(1, 2) + m_turb_cool; % Kg / s
PMAT(6, 3) = A_circ(0, 102.26E-3); % m^2 (NPS 4 Sch 40 Pipe)
PMAT(6, 4) = fzero(@ (x) (integral(Cp_H2, PMAT(2, 4), PMAT(5, 4), "AbsTol", 1E-15) - (PWR_Turb ./ PMAT(6, 2)) - integral(Cp_H2, PMAT(2, 4), x, "AbsTol", 1E-15)), [PMAT(2, 4) PMAT(5, 4)]); % K


% Setup Point 7 - Nozzle Outlet
PMAT(7, 1) = P0 / (T0_by_T(M_e, gamma_H2)^(gamma_H2 / (gamma_H2 - 1))); % Pa
PMAT(7, 2) = (1 - M_y) * PMAT(4, 2); % Kg / s
PMAT(7, 3) = Ae; % m^2 (NPS 4 Sch 40 Pipe)
PMAT(7, 4) = T0 / T0_by_T(M_e, gamma_H2); % K

% Calculate Rocket Performance
Thrust = PMAT(7, 2) * M_e * sqrt(gamma_H2 * R_H2 * PMAT(7, 4)) + Ae * PMAT(7, 1); % N
ISP_vac = Thrust / (PMAT(7, 2) * 9.81); % s

% Calculate Thermodynamic Efficiencies
eta_th = 1 - (PMAT(6, 4) ./ PMAT(1, 4));
eta_ac =  PWR_net / Q_c;

% Output Calculated Values
% Input Values
Q_c
P_init
T_init
m_dot
M_y
Ae_by_A_t
printf("\n");

% Temperature Limit Input Data
T_mar
Tmax_Al2O3
Tmax_MgO
Tmax_UC
Tmax_turbin
printf("\n");

% Chamber and Nozzle Data
m_UC
t_PV
R_in
R_t
R_e
Z_C
Z_D
Z_sep
Thrust
ISP_vac
printf("\n");

% Reactor Temperature Data
TmaxM
TmaxNorm
TmaxZ
printf("\n");

% Turbomachinery Data
m_turb_cool
PWR_Turb
PWR_Pmp
PWR_net
printf("\n");

% Thermal Efficiency Data
eta_th
eta_ac
printf("\n");

% Reactor Pressure - Mass Flow - Area - Temperature Data
PMAT

