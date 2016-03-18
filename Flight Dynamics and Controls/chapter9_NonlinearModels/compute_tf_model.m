function [T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,P)
% x_trim is the trimmed state,
% u_trim is the trimmed input

% Constants from chapter 3
alpha = atan(P.w0/P.u0);

Gamma = P.Jx*P.Jz - P.Jxz^2;
Gamma_1 = (P.Jxz*(P.Jx- P.Jy + P.Jz))/Gamma;
Gamma_2 = (P.Jz*(P.Jz-P.Jy) + P.Jxz^2)/Gamma;
Gamma_3 = P.Jz/Gamma;
Gamma_4 = P.Jxz/Gamma;
Gamma_5 = (P.Jz-P.Jx)/P.Jy;
Gamma_6 = P.Jxz/P.Jy;
Gamma_7 = ((P.Jx-P.Jy)*P.Jx + P.Jxz^2)/Gamma;
Gamma_8 = P.Jx/Gamma;

% Constants from chapter 5
Cp_0        = Gamma_3*P.C_ell_0       + Gamma_4*P.C_n_0;
Cp_beta     = Gamma_3*P.C_ell_beta    + Gamma_4*P.C_n_beta;
Cp_p        = Gamma_3*P.C_ell_p       + Gamma_4*P.C_n_p;
Cp_r        = Gamma_3*P.C_ell_r       + Gamma_4*P.C_n_r;
Cp_delta_a  = Gamma_3*P.C_ell_delta_a + Gamma_4*P.C_n_delta_a;
Cp_delta_r  = Gamma_3*P.C_ell_delta_r + Gamma_4*P.C_n_delta_r;
Cr_0        = Gamma_4*P.C_ell_0       + Gamma_8*P.C_n_0;
Cr_beta     = Gamma_4*P.C_ell_beta    + Gamma_8*P.C_n_beta;
Cr_p        = Gamma_4*P.C_ell_p       + Gamma_8*P.C_n_p;
Cr_r        = Gamma_4*P.C_ell_r       + Gamma_8*P.C_n_r;
Cr_delta_a  = Gamma_4*P.C_ell_delta_a + Gamma_8*P.C_n_delta_a;
Cr_delta_r  = Gamma_4*P.C_ell_delta_r + Gamma_8*P.C_n_delta_r;

% T_phi_delta_a constants
a_phi1  = -(1/2)*P.rho*(P.Va0^2)*P.S_wing*P.b*Cp_p*(P.b/(2*P.Va0));
a_phi2  =  (1/2)*P.rho*(P.Va0^2)*P.S_wing*P.b*Cp_delta_a;

% T_chi_phi constants
Va_trim = P.Va0*cos(P.gamma0); %ground speed with no wind

% T_theta_delta_e constants
a_theta1 = -(P.rho*P.Va0^2*P.c*P.S_wing/(2*P.Jy))*P.C_m_q*(P.c/(2*P.Va0));
a_theta2 = -(P.rho*P.Va0^2*P.c*P.S_wing/(2*P.Jy))*P.C_m_alpha;
a_theta3 = -(P.rho*P.Va0^2*P.c*P.S_wing/(2*P.Jy))*P.C_m_delta_e;

% T_h_Va constants
theta_trim = P.theta0;

% T_Va_delta_t constants
a_V1 = P.rho*Va_trim*P.S_wing/(P.mass)*(P.C_D_0 + P.C_D_alpha*alpha + P.C_D_delta_e*u_trim(1))...
        + (P.rho*P.S_prop/P.mass)*P.C_prop*Va_trim;
a_V2 = (P.rho*P.S_prop/P.mass)*P.C_prop*(P.k_motor^2)*u_trim(4);
a_V3 = P.gravity*cos(P.theta0 - alpha);

% T_v_delta_r
a_beta1 = -(P.rho*P.Va0*P.S_wing/(2*P.M))*P.C_Y_beta;
a_beta2 =  (P.rho*P.Va0*P.S_wing/(2*P.M))*P.C_Y_delta_r;

% define transfer functions
T_phi_delta_a   = tf([a_phi2],[1,a_phi1,0]);
T_chi_phi       = tf([P.gravity/Va_trim],[1,0]);
T_theta_delta_e = tf(a_theta3,[1,a_theta1,a_theta2]);
T_h_theta       = tf([Va_trim],[1,0]);
T_h_Va          = tf([theta_trim],[1,0]);
T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
T_v_delta_r     = tf([Va_trim*a_beta2],[1,a_beta1]);

