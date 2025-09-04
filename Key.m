% Author: Khushi Piparava
% Purpose:
%   Exhaustive search for a double-wedge supersonic airfoil that maximizes
%   L/D while meeting minimum lift and pitch-moment constraints at M = 2.6.
%   Geometry is defined by two linear panels per side with vertex locations
%   (xu, tu) on the upper surface and (xl, tl) on the lower surface.
%
% Design spec (matches assignment brief):
%   - Design Mach number: M1 = 2.6
%   - Angle of attack: a = 5 deg
%   - Objective: maximize Cl/Cd
%   - Constraints:
%       * tu + tl = 0.1 (10 percent thickness)
%       * |Cm_LE| <= 0.05
%       * Cl >= 0.2
% Reference: shock-expansion method for 2D supersonic flow
% Note: This script expects the following helper functions on the MATLAB path:
%   oswbeta, nswr_newmach, nswr_pressure, pm, stag2stat
%
% Performance note:
%   The full 3D grid can be very large. Start with coarser steps, then refine.

format short
clc

% -----------------------------
% 1) Setup and containers
% -----------------------------
chord_length = 1;
M1 = 2.6;
a  = 5;        % deg
g  = 1.4;

% Search grids
tu_grid = 0:0.002:0.1;   % upper thickness
xu_grid = 0:0.002:1;     % upper vertex x/c
xl_grid = 0:0.002:1;     % lower vertex x/c

% Results storage
Cl_values = [];
Cd_values = [];
Cm_values = [];
tu_values = [];
tl_values = [];
xu_values = [];
xl_values = [];

% -----------------------------
% 2) Exhaustive search
% -----------------------------
for tu = tu_grid
    tl = 0.1 - tu;                 % enforce fixed total thickness
    if tl < 0, continue; end       % guard against roundoff

    for xu = xu_grid
        % Geometry angles in degrees
        theta1 = atand( max(tu,eps) / max(xu,eps) );        % upper LE panel
        theta3 = atand( max(tu,eps) / max(1 - xu,eps) );    % upper TE panel

        for xl = xl_grid
            theta2 = atand( max(tl,eps) / max(xl,eps) );        % lower LE panel
            theta4 = atand( max(tl,eps) / max(1 - xl,eps) );    % lower TE panel

            % -----------------------------
            % 2a) Lower surface aerodynamics
            %     LE: oblique shock, vertex: expansion
            % -----------------------------
            beta_l = oswbeta(M1, theta1 + theta2, g);               % shock angle [deg]
            M_n1   = M1 * sin(deg2rad(beta_l));
            M_n2   = nswr_newmach(M_n1);
            M_l2   = M_n2 / sin(deg2rad(beta_l - (theta1 + theta2)));
            pl2_p1 = nswr_pressure(M_n1);
            Cpl1   = (pl2_p1 - 1) / (0.5 * g * M1^2);               % panel 1

            M_l3   = pm(M_l2, theta2 + theta4, g);                  % expansion to panel 2
            P0_pl2 = stag2stat(M_l2);
            P0_pl3 = stag2stat(M_l3);
            pl3_pl2 = P0_pl2 / P0_pl3;
            pl3_p1  = pl3_pl2 * pl2_p1;
            Cpl2    = (pl3_p1 - 1) / (0.5 * g * M1^2);              % panel 2

            % -----------------------------
            % 2b) Upper surface aerodynamics
            %     If a > theta1, LE is expansion. Else shock. Vertex is expansion.
            % -----------------------------
            if a > theta1
                % Expansion at upper LE
                M_u2   = pm(M1, a - theta1, g);
                P0_pu1 = stag2stat(M1);
                P0_pu2 = stag2stat(M_u2);
                pu2_pu1 = P0_pu1 / P0_pu2;
                Cpu1    = (pu2_pu1 - 1) / (0.5 * g * M1^2);
            elseif a == theta1
                M_u2 = M1;
                Cpu1 = 0;
                pu2_pu1 = 1;
            else
                % Oblique shock at upper LE
                beta_u  = oswbeta(M1, theta1, g);
                M_n1_u  = M1 * sin(deg2rad(beta_u));
                M_n2_u  = nswr_newmach(M_n1_u);
                M_u2    = M_n2_u / sin(deg2rad(beta_u - theta1));
                pu2_pu1 = nswr_pressure(M_n1_u);
                Cpu1    = (pu2_pu1 - 1) / (0.5 * g * M1^2);
            end

            % Expansion at upper vertex
            M_u3   = pm(M_u2, theta2 + theta4, g);
            P0_pu2 = stag2stat(M_u2);
            P0_pu3 = stag2stat(M_u3);
            pu3_pu2 = P0_pu2 / P0_pu3;
            pu3_p1  = pu3_pu2 * pu2_pu1;
            Cpu2    = (pu3_p1 - 1) / (0.5 * g * M1^2);

            % -----------------------------
            % 2c) Force and moment coefficients
            %     Sign conventions follow panel orientations
            % -----------------------------
            Cn = Cpl1*xl + Cpl2*(1 - xl) - Cpu1*xu - Cpu2*(1 - xu);
            Ca =  Cpu1*xu*tan(deg2rad(theta1)) ...
                + Cpu2*(1 - xu)*tan(-deg2rad(theta3)) ...
                - Cpl1*xl*tan(-deg2rad(theta2)) ...
                - Cpl2*(1 - xl)*tan( deg2rad(theta4));

            Cl = Cn*cosd(a) - Ca*sind(a);
            Cd = Cn*sind(a) + Ca*cosd(a);

            % Leading-edge moment coefficient about x = 0
            Cm =  Cpu1*(xu^2)/2 ...
                + Cpu2*((1 - xu^2))/2 ...
                - Cpl1*(xl^2)/2 ...
                - Cpl2*((1 - xl^2))/2 ...
                + Cpu1*(xu*tan(deg2rad(theta1)))^2/2 ...
                + Cpu2*((1 - xu)*tan(-deg2rad(theta3)))^2/2 ...
                - Cpl1*(xl*tan(-deg2rad(theta2)))^2/2 ...
                - Cpl2*((1 - xl)*tan( deg2rad(theta4)))^2/2;

            % -----------------------------
            % 2d) Constraint check and save
            % -----------------------------
            if (Cl >= 0.2) && (abs(Cm) <= 0.05)
                Cl_values(end+1) = real(Cl); %#ok<*AGROW>
                Cd_values(end+1) = real(Cd);
                Cm_values(end+1) = real(Cm);
                tu_values(end+1) = tu;
                tl_values(end+1) = tl;
                xu_values(end+1) = xu;
                xl_values(end+1) = xl;
            end
        end
    end
end

% -----------------------------
% 3) Pick best design by max L/D
% -----------------------------
if isempty(Cl_values)
    error('No feasible designs found. Widen your grid or relax steps.');
end

[~, max_index] = max(Cl_values ./ Cd_values);
optimal_xu = xu_values(max_index);
optimal_xl = xl_values(max_index);
optimal_tu = tu_values(max_index);
optimal_tl = tl_values(max_index);
optimal_Cl = Cl_values(max_index);
optimal_Cd = Cd_values(max_index);
optimal_Cm = Cm_values(max_index);

fprintf('Best design:\n xu=%.4f, xl=%.4f, tu=%.4f, tl=%.4f\n', ...
    optimal_xu, optimal_xl, optimal_tu, optimal_tl);
fprintf(' Cl=%.4f, Cd=%.4f, Cm=%.4f, L/D=%.2f\n', ...
    optimal_Cl, optimal_Cd, optimal_Cm, optimal_Cl/optimal_Cd);

% -----------------------------
% 4) Plot the best airfoil
% -----------------------------
X = [0, optimal_xu, 1, optimal_xl, 0];
Y = [0, optimal_tu, 0, -optimal_tl, 0];

figure;
plot(X, Y, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Chord Length');
ylabel('Thickness');
title('Optimal Airfoil Shape');
xlim([0, chord_length]);
ylim([-0.4, 0.4]);
grid on;
