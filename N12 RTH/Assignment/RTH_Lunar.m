% This Octave/MATLAB Script Performs Thermohydraulic Calculations
% For a NERVA Type NTP on Luna

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

Q_c_max = 2.25E9; % W
Q_ele_max = 0.65E6; % W
Q_c = 1.82E9; % W
Num_ele = Q_c_max / Q_ele_max;
Q_ele = Q_c / Num_ele; % W
H_ele = 0.89; % m
H_ex = 0.93463; % m

T_mar = 40; % K
SF_PV = 1.5;
Sheath_W = 1; % Units of t_PV

P_init = 4.5053E6; % Pa
T_init = 1373; % K
m_dot = 87.25; % Kg / s

% Define Original Space Reactor Geometry
C_deg = deg2rad(35); % rad
D_deg = deg2rad(15); % rad
Ae_by_A_t = 100;

m_UC = 1.3476e+04; % Kg
t_PV = 0.020000; % m
R_in = 0.6299; % m
R_t = 0.1020; % m
R_e = 1.0204; % m
Z_C = 0.7538; % m
Z_D = 3.4275; % m
Z_sep = 1.1441e-04; % m
A_rad = (pi .* (R_t + R_e) .* Z_C ./ cos(D_deg)); % m^2

% Define Turbine Data
Tmax_turbin = 1373; % K (Maximum Uncooled Turbine Blade Limit)
Num_Comp_Stages = 3;
Num_Turb_per_Module = 2;
PRATIO_Turb = (1.61^3)^Num_Turb_per_Module; % Raised to the Number of Stages, Raised to Number of Turbines
Tmin = 20; % K (Minimum Fluid Loop Temperature)

% Define H2 Material Data Functions
% Load H2 Material Data
Cp_H2_dat = csvread("H2_Cp_5MPa.csv");
mu_H2_dat = csvread("H2_mu_5MPa.csv");
k_H2_dat = csvread("H2_k_5MPa.csv");
Pr_H2_dat = csvread("H2_Pr_5MPa.csv");

% Isobaric Specific Heat vs Temperature
Cp_H2 = @ (T) interp1(Cp_H2_dat(:, 1), Cp_H2_dat(:, 2), T, "cubic", "extrap"); % J / Kg * K

% Dynamic Viscosity vs Temperature
mu_H2 = @ (T) interp1(mu_H2_dat(:, 1), mu_H2_dat(:, 2), T, "cubic", "extrap"); % Pa * s

% Thermal Conductivity vs Temperature
k_H2 = @ (T) interp1(k_H2_dat(:, 1), k_H2_dat(:, 2), T, "cubic", "extrap"); % W / m * K

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

% Black Body Emissivity
e_Al2O3 = 0.8;

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
% 1 - Chamber Outlet / Turbine Inlet
% 2 - Turbine Outlet / Pump Inlet
% 3 - Pump Outlet / Sheath Entry
% 4 - Sheath Exit / Reactor Plenum
% 5 - Reactor Outlet / Chamber Inlet
PMAT = zeros(5, 4);


% Setup Point 1 - Chamber Outlet / Turbine Inlet
PMAT(1, 1) = P_init; % Pa
PMAT(1, 2) = m_dot; % Kg / s
PMAT(1, 3) = A_circ(0, 102.26E-3); % m^2 (NPS 4 Sch 40 Pipe)
PMAT(1, 4) = T_init; % K

% Calculate Turbine Output
PWR_Turb = 0.9 .* 0.7171 .* PMAT(1, 2) .* Cp_H2(PMAT(1, 4)) .* PMAT(1, 4) .* ((PRATIO_Turb.^((gamma_H2 - 1) ./ gamma_H2)) - 1); % W


% Setup Point 2 - Turbine Outlet / Pump Inlet
PMAT(2, 1) = PMAT(1, 1) / PRATIO_Turb; % Pa
PMAT(2, 2) = PMAT(1, 2); % Kg / s
PMAT(2, 3) = A_circ(R_t + t_PV, R_t + (1 + Sheath_W) * t_PV); % m^2
PMAT(2, 4) = fzero(@ (x) (integral(Cp_H2, Tmin, PMAT(1, 4), "AbsTol", 1E-15) - (PWR_Turb ./ PMAT(2, 2)) - integral(Cp_H2, Tmin, x, "AbsTol", 1E-15)), [Tmin PMAT(1, 4)]); % K

##% Calculate Radiator Data Using Thin Wall Approximation
##f_rad = @ (x) (Cp_H2(x) ./ (x.^4 - Tamb.^4));
##Fv = 1 + sin(D_deg); % Radiative View Factor
##sigma_sb = 5.670374E-8; % W / m^2 * K^4 (Stefan-Boltzmann Constant)
##T_rad_out = fzero(@ (x) (integral(f_rad, x, PMAT(3, 4), "AbsTol", 1E-15) - (Fv .* e_Al2O3 .* sigma_sb .* A_rad ./ PMAT(3, 2))), [Tmin PMAT(3, 4)]); % K
##DP_rad = fzero(@ (x) (x - (((PMAT(3, 2) ./ PMAT(3, 3)).^2) .* (rho(PMAT(3, 1) - x, R_H2, T_rad_out).^-1 - rho(PMAT(3, 1), R_H2, PMAT(3, 4)).^-1))), [0 1E5]); % Pa

% Calculate Pump Input Work at Pump Outlet Pressure of 5 MPa and Compression Ratio
PWR_Pmp = (PMAT(2, 2) .* Cp_H2(PMAT(2, 4)) .* PMAT(2, 4) .* (((5E6 ./ PMAT(2, 1)).^((gamma_H2 - 1) ./ gamma_H2)) - 1)) / (0.9 * 0.7171); % W
PRATIO_Pmp = power((5E6 ./ PMAT(2, 1)), Num_Comp_Stages.^-1);


% Setup Point 3 - Pump Outlet / Sheath Inlet
PMAT(3, 1) = 5E6; % Pa
PMAT(3, 2) = PMAT(2, 2); % Kg / s
PMAT(3, 3) = A_circ(R_t + t_PV, R_t + (1 + Sheath_W) * t_PV); % m^2
PMAT(3, 4) = fzero(@ (x) (integral(Cp_H2, Tmin, PMAT(2, 4), "AbsTol", 1E-15) + (PWR_Pmp ./ PMAT(2, 2)) - integral(Cp_H2, Tmin, x, "AbsTol", 1E-15)), [Tmin Tmax_turbin]); % K

% Calculate Net Eletrical Work and Cycle Efficiency
PWR_net = PWR_Turb - PWR_Pmp; % W
eta_th = 1 - (PMAT(3, 4) ./ PMAT(1, 4));
PWR_eta = Q_c * eta_th; % W
eta_ac = PWR_net / Q_c;


% Calculate Regeneratively Cooled Nozzle Data Using Conduction Approximation
Q_rgen_cond = 15 * pi * (R_in + R_t) * Z_C / cos(C_deg) * PMAT(1, 4) / t_PV; % W (A-286 Conductivity is 15 W / m * K)
T_rgen_cond = fzero(@ (x) (integral(Cp_H2, PMAT(3, 4), x, "AbsTol", 1E-15) - (Q_rgen_cond / PMAT(3, 2))), [PMAT(3, 4) Tmax_turbin]); % K

% Setup and Calculate Maximum Pressure Loss in Nozzle Sheath
##DP_rgen = 0;
DP_rgen = fzero(@ (x) (x - (((PMAT(3, 2) ./ PMAT(3, 3)).^2) .* (rho(PMAT(3, 1), R_H2, T_rgen_cond).^-1 - rho(PMAT(3, 1), R_H2, PMAT(3, 4)).^-1))), [0 1E6]); % Pa

% Setup Point 4 - Sheath Outlet / Reactor Plenum
PMAT(4, 1) = PMAT(3, 1) - DP_rgen; % Pa
PMAT(4, 2) = PMAT(3, 2); % Kg / s
PMAT(4, 3) = A_circ(R_in + t_PV, R_in + (1 + Sheath_W) * t_PV); % m^2
PMAT(4, 4) = T_rgen_cond; % K


% Setup Reactor Core Calculation
% Initialize Reactor Temperature and Power Profile Container
div = 1024;
TBCSM = zeros(div, 4);
DP_c = 0; % Pa

% Calculate Reactor Element Channel Data and Setup Helper Functions
A_chnl = A_circ(0, D_chnl / 2); % m^2
A_flow = 19 * A_chnl; % m^2
M_chnl = (PMAT(4, 2) * Q_ele) / (19 * Q_c); % Kg / s
Q_fc = @ (z) (H_ele ./ div) .* ((Q_ele ./ 19) .* cos(pi .* ((H_ele / 2) - z) ./ H_ex)); % W
HTC  = @ (T) ((k_H2(T) ./ D_chnl) .* 0.205 .* (f_Re(M_chnl, 0, D_chnl / 2, T).^0.8) .* (Pr_H2(T).^0.4)); % W / m^2 * K

% Calculate Inlet Initial Temperatures
TBCSM(1, 1) = PMAT(4, 4); % K
TBCSM(1, 2) = fzero(@ (x) (TBCSM(1, 1) - x + (Q_fc(H_ele ./ div) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(1, 1)) .* ((x ./ TBCSM(1, 1)).^-0.55))))), [PMAT(4, 4) Tmax_MgO]); % K
##TBCSM(1, 2) = (TBCSM(1, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(1, 1))))); % K
TBCSM(1, 3) = ((Q_fc(H_ele ./ div) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(1, 2)))) + TBCSM(1, 2); % K
TBCSM(1, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(1, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(1, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(1, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(1, 3)); % K
TBCSM(2, 1) = TBCSM(1, 1) + Q_fc(H_ele ./ div) ./ (M_chnl .* Cp_H2(TBCSM(1, 1))); % K
DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(1, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(1, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa

% Calculate Full TBCSM Matrix
for i = 2:(div - 1)
  TBCSM(i, 2) = fzero(@ (x) (TBCSM(i, 1) - x + (Q_fc(i .* H_ele ./ div) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(i, 1)) .* ((x ./ TBCSM(i, 1)).^-0.55))))), [PMAT(4, 4) Tmax_MgO]); % K
##  TBCSM(i, 2) = (TBCSM(i, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(i, 1))))); % K
  TBCSM(i, 3) = ((Q_fc(i .* H_ele ./ div) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(i, 2)))) + TBCSM(i, 2); % K
  TBCSM(i, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(i, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(i, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(i, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(i, 3)); % K
  TBCSM(i + 1, 1) = TBCSM(i, 1) + Q_fc(i .* H_ele ./ div) ./ (M_chnl .* Cp_H2(TBCSM(i, 1))); % K
  DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(i, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(i, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa
endfor

% Calculate Core Outlet Values for TBCSM Matrix
TBCSM(div, 2) = fzero(@ (x) (TBCSM(div, 1) - x + (Q_fc(H_ele) ./ ((pi .* D_chnl .* H_ele ./ div) .* (HTC(TBCSM(div, 1)) .* ((x ./ TBCSM(div, 1)).^-0.55))))), [PMAT(4, 4) Tmax_MgO]); % K
##TBCSM(div, 2) = (TBCSM(div, 1) + (Q_fc(H_ele ./ div) ./ (pi .* D_chnl .* H_ele ./ div .* HTC(TBCSM(div, 1))))); % K
TBCSM(div, 3) = ((Q_fc(H_ele) .* log(D_hole ./ D_chnl)) ./ (2 .* pi .* H_ele ./ div .* k_MgO(TBCSM(div, 2)))) + TBCSM(div, 2); % K
TBCSM(div, 4) = ((Q_fc(H_ele ./ div) ./ (H_ele ./ div)) .* ((1 ./ (4 .* pi .* k_UC(TBCSM(div, 3)))) + log(D_sep ./ D_hole) .* ((1 ./ (2 .* pi .* k_MgO(TBCSM(div, 3)))) + (((0.5 .* D_hole).^2) ./ (2 .* k_UC(TBCSM(div, 3)) .* A_circ(0.5 .* D_hole, 0.5 .* D_sep))))) + TBCSM(div, 3)); % K
DP_c += ((H_ele ./ div) .* (f_dw(f_Re(M_chnl, 0, D_chnl / 2, TBCSM(div, 1)) .* ((M_chnl ./ A_chnl).^2)) ./ (8 .* rho((PMAT(4, 1) - DP_c), R_H2, TBCSM(div, 1)) .* (A_chnl ./ (pi .* D_chnl))))); % Pa

% Find Maximum Temperatures Locations
[TmaxM, TmaxInd] = max(TBCSM);
TmaxNorm = TmaxM ./ [Tmax_turbin Tmax_MgO Tmax_UC Tmax_UC];
TmaxZ = (H_ele / 2) .- TmaxInd .* (H_ele / div); % m

% Calculate Core Pressure Loss
DP_c += fzero(@ (x) (x - (((M_chnl ./ A_chnl).^2) .* (rho(PMAT(4, 1) - x, R_H2, TBCSM(1024, 1)).^-1 - rho(PMAT(4, 1), R_H2, PMAT(4, 4)).^-1))), [0 1.75E6]); % Pa


% Setup Point 5 - Reactor Outlet / Chamber Inlet
PMAT(5, 1) = PMAT(4, 1) - DP_c; % Pa
PMAT(5, 2) = PMAT(4, 2); % Kg / s
PMAT(5, 3) = A_circ(0, R_in); % m^2 (NPS 4 Sch 40 Pipe)
PMAT(5, 4) = TBCSM(div, 1); % K

% Output Calculated Values
% Input Values
Q_c
P_init
T_init
m_dot
printf("\n");

% Temperature Limit Input Data
T_mar
Tmax_Al2O3
Tmax_MgO
Tmax_UC
Tmax_turbin
printf("\n");

% Chamber and Nozzle Data
t_PV
R_in
R_t
R_e
Z_C
Z_D
Z_sep
printf("\n");

% Reactor Temperature Data
TmaxM
TmaxNorm
TmaxZ
printf("\n");

% Turbomachinery Module Data
PRATIO_Pmp
Num_Comp_Stages
PRATIO_Turb
Num_Turb_per_Module
printf("\n");

% Plant Electrical Data
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

