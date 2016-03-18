function [T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,y_trim,P)
% x_trim is the trimmed state,
% u_trim is the trimmed input

Va_trim = x_trim(4);
theta_trim = x_trim(8);

alpha_trim = y_trim(2);
gamma_trim = y_trim(3);

delta_e_trim = u_trim(1);
delta_t_trim = u_trim(4);

% define intermediate variables
G = P.Jx*P.Jz-(P.Jxz)^2;
G1 = P.Jxz*(P.Jx-P.Jy+P.Jz) / G;
G2 = (P.Jz*(P.Jz-P.Jy)+(P.Jxz^2)) / G;
G3 = P.Jz / G;
G4 = P.Jxz / G;
G5 = (P.Jz-P.Jx)/P.Jy;
G6 = P.Jxz/P.Jy;
G7 = ((P.Jx-P.Jy)*P.Jx+(P.Jxz)^2)/G;
G8 = P.Jx / G;
ms = P.mass;
g = P.gravity;
CL0 = P.C_L_0;
CLa = P.C_L_alpha;
CLq = P.C_L_q;
CLde = P.C_L_delta_e;
Cd0 = P.C_D_0;
Cda = P.C_D_alpha;
Cdp = P.C_D_p;
Cdq = P.C_D_q;
Cdde = P.C_D_delta_e;
Cy0 = P.C_Y_0;
Cyb = P.C_Y_beta;
Cyp = P.C_Y_p;
Cyr = P.C_Y_r;
Cyda = P.C_Y_delta_a;
Cydr = P.C_Y_delta_r;
Cl0 = P.C_ell_0;
Clb = P.C_ell_beta;
Clp = P.C_ell_p;
Clr = P.C_ell_r;
Clda = P.C_ell_delta_a;
Cldr = P.C_ell_delta_r;
Cm0 = P.C_m_0;
Cma = P.C_m_alpha;
Cmq = P.C_m_q;
Cmde = P.C_m_delta_e;
Cn0 = P.C_n_0;
Cnb = P.C_n_beta;
Cnp = P.C_n_p;
Cnr = P.C_n_r;
Cnda = P.C_n_delta_a;
Cndr = P.C_n_delta_r;
komega = P.k_Omega;
ktp = P.k_T_P;
AR = (P.b)^2/P.S_wing;
Cd = Cdp + (CL0+CLa*alpha_trim)^2/(pi()*P.e*AR);
%{
sigma = (1+exp(-P.M*(alpha-P.alpha0))+exp(P.M*(alpha+P.alpha0))) / ...
            ((1+exp(-P.M*(alpha-P.alpha0)))*(1+exp(P.M*(alpha+P.alpha0))));
CL = (1-sigma)*(CL0+CLa*alpha)+sigma*(2*sign(alpha)*(sin(alpha))^2*cos(alpha));
Cx = -Cd*cos(alpha)+CL*sin(alpha);
Cxq = -Cdq*cos(alpha)+CLq*sin(alpha);
Cxde = -Cdde*cos(alpha)+CLde*sin(alpha);
Cz = -Cd*sin(alpha)-CL*cos(alpha);
Czq = -Cdq*sin(alpha)-CLq*cos(alpha);
Czde = -Cdde*sin(alpha)-CLde*cos(alpha);
%}

% begin tf model coefficients
Cp0 = G3*Cl0+G4*Cn0;
Cpb = G3*Clb+G4*Cnb;
Cpp = G3*Clp+G4*Cnp;
Cpr = G3*Clr+G4*Cnr;
Cpda = G3*Clda+G4*Cnda;
Cpdr = G3*Cldr+G4*Cndr;
Cr0 = G4*Cl0+G8*Cn0;
Crb = G4*Clb+G8*Cnb;
Crp = G4*Clp+G8*Cnp;
Crr = G4*Clr+G8*Cnr;
Crda = G4*Clda+G8*Cnda;
Crdr = G4*Cldr+G8*Cndr;
Va = P.Va0;

a_phi1 = -P.rho * Va^2 * P.S_wing * P.b * Cpp * P.b / (2*2*Va);
a_phi2 = P.rho * Va^2 * P.S_wing * P.b * Cpda / 2;
a_theta1 = - (P.rho * Va^2 * P.c * P.S_wing) / (2*P.Jy) * Cmq * P.c / (2 * Va);
a_theta2 = - (P.rho * Va^2 * P.c * P.S_wing) / (2*P.Jy) * Cma;
a_theta3 = (P.rho * Va^2 * P.c * P.S_wing) / (2*P.Jy) * Cmde;
a_V1 = P.rho * Va_trim * P.S_wing / (P.mass) * (Cd0 + Cda * alpha_trim * Cdde * delta_e_trim) + P.rho * P.S_prop / P.mass * P.C_prop * Va_trim;
a_V2 = P.rho * P.S_prop / P.mass * P.C_prop * P.k_motor^2 * delta_t_trim;
a_V3 = P.gravity * cos(theta_trim - gamma_trim);
a_beta1 = - P.rho * Va * P.S_wing / (2 * P.mass) * Cyb;
a_beta2 = P.rho * Va * P.S_wing / (2 * P.mass) * Cydr;

    
% define transfer functions
T_phi_delta_a   = tf([a_phi2],[1,a_phi1,0]);
T_chi_phi       = tf([P.gravity/Va_trim],[1,0]);
T_theta_delta_e = tf(a_theta3,[1,a_theta1,a_theta2]);
T_h_theta       = tf([Va_trim],[1,0]);
T_h_Va          = tf([theta_trim],[1,0]);
T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
T_v_delta_r     = tf([Va_trim*a_beta2],[1,a_beta1]);

