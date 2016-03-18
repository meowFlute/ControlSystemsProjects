function [T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,y_trim,P)
% x_trim is the trimmed state,
% u_trim is the trimmed input

% add stuff here

    jx = P.Jx;
    jy = P.Jy;
    jz = P.Jz;
    jxz= P.Jxz;
    
    gamma = jx*jz-jxz^2;
    gamma1 = (jxz*(jx-jy+jz))/gamma;
    gamma2 = (jz*(jz-jy)+jxz^2)/gamma;
    gamma3 = jz/gamma;
    gamma4 = jxz/gamma;
    gamma5 = (jz-jx)/jy;
    gamma6 = jxz/jy;
    gamma7 = ((jx-jy)*jx+jxz^2)/gamma;
    gamma8 = jx/gamma;
    
    rho = P.rho;
    Va = P.Va0;
    Va_trim = P.Va0;
    theta_trim = P.theta0;
    gamma_trim = y_trim(3);
    
    c = P.c;
    S = P.S_wing;
    b = P.b;
    Cp0 = gamma3*P.C_ell_0+gamma4*P.C_n_0;
    Cpb = gamma3*P.C_ell_beta+gamma4*P.C_n_beta;
    Cpp = gamma3*P.C_ell_p+gamma4*P.C_n_p;
    Cpr = gamma3*P.C_ell_r+gamma4*P.C_n_r;
    Cpda = gamma3*P.C_ell_delta_a+gamma4*P.C_n_delta_a;
    Cpdr = gamma3*P.C_ell_delta_r+gamma4*P.C_n_delta_r;
    
    Cr0 = gamma4*P.C_ell_0+gamma8*P.C_n_0;
    Crb = gamma4*P.C_ell_beta+gamma8*P.C_n_beta;
    Crp = gamma4*P.C_ell_p+gamma8*P.C_n_p;
    Crr = gamma4*P.C_ell_r+gamma8*P.C_n_r;
    Crda = gamma4*P.C_ell_delta_a+gamma8*P.C_n_delta_a;
    Crdr = gamma4*P.C_ell_delta_r+gamma8*P.C_n_delta_r;
    
    a_phi1 = -1/2*rho*Va^2*S*b*Cpp*b/Va/2;
    a_phi2 = 1/2*rho*Va^2*S*b*Cpda;
    
    a_beta1 = -(rho*Va*S)/(2*P.mass)*P.C_Y_beta;
    a_beta2 = (rho*Va*S)/(2*P.mass)*P.C_Y_delta_r;
    
    a_theta1 = -(rho*Va^2*c*S)/(2*P.Jy)*P.C_m_q*c/(2*Va);
    a_theta2 = -(rho*Va^2*c*S)/(2*P.Jy)*P.C_m_alpha;
    a_theta3 = -(rho*Va^2*c*S)/(2*P.Jy)*P.C_m_delta_e;
    
    a_V1 = (rho*Va_trim*S)/P.mass*(P.C_D_0+P.C_D_alpha*y_trim(2)+P.C_D_delta_e*u_trim(1))+(rho*P.S_prop)/P.mass*P.C_prop*Va_trim;
    a_V2 = (rho*P.S_prop)/P.mass*P.C_prop*P.k_motor^2*u_trim(4);
    a_V3 = P.g*cos(theta_trim-gamma_trim);
    
% define transfer functions
T_phi_delta_a   = tf([a_phi2],[1,a_phi1,0]);
T_chi_phi       = tf([P.gravity/Va_trim],[1,0]);
T_theta_delta_e = tf(a_theta3,[1,a_theta1,a_theta2]);
T_h_theta       = tf([Va_trim],[1,0]);
T_h_Va          = tf([theta_trim],[1,0]);
T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
T_v_delta_r     = tf([Va_trim*a_beta2],[1,a_beta1]);

