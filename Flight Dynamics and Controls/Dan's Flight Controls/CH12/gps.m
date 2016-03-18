% gps.m
%   Compute the output of gps sensor
%
%  Revised:
%   3/5/2010 - RB 
%   5/14/2010 - RB

function y = gps(uu, P)

    % relabel the inputs
    Va      = uu(1);
%    alpha   = uu(2);
%    beta    = uu(3);
    wn      = uu(4);
    we      = uu(5);
%    wd      = uu(6);
    pn      = uu(7);
    pe      = uu(8);
    pd      = uu(9);
%    u       = uu(10);
%    v       = uu(11);
%    w       = uu(12);
%    phi     = uu(13);
%    theta   = uu(14);
    psi     = uu(15);
%    p       = uu(16);
%    q       = uu(17);
%    r       = uu(18);
    t       = uu(19);
    
    persistent vn;
    persistent ve;
    persistent vh;
    
    persistent pos_north_d1;
    persistent pos_east_d1;
    
    if t == 0
        vn = 0;
        ve = 0;
        vh = 0;
        pos_north_d1 = 0;
        pos_east_d1 = 0;
    end
          
    % construct North, East, and altitude GPS measurements
    y_gps_n = pn + vn;
    y_gps_e = pe + ve; 
    y_gps_h = pd + vh; 
    
    % construct groundspeed and course measurements
    V_north = (y_gps_n - pos_north_d1) / P.Ts;
    V_east = (y_gps_e - pos_east_d1) / P.Ts;
    V_g = sqrt(V_north^2 + V_east^2);
    sigma_V_g = sqrt((V_north^2*P.sigma_n_gps^2 + V_east^2*P.sigma_e_gps^2));
    if sigma_V_g ~= 0
        sigma_V_g = sigma_V_g / sqrt((V_north^2+V_east^2)^2);
    end
    %sigma_V_chi = sqrt((V_north^2*P.sigma_n_gps^2 + V_east^2*P.sigma_e_gps^2) / (V_north^2+V_east^2)^2);
    sigma_V_chi = sigma_V_g;
    if sigma_V_chi ~= 0
        sigma_V_chi = sigma_V_chi / V_g;
    end
    
    y_gps_Vg     = sqrt((Va * cos(psi) + wn)^2 + (Va * sin(psi) + we)^2) + sigma_V_g * randn();
    y_gps_course = atan2(Va * sin(psi) + we, Va * cos(psi) + wn) + sigma_V_chi * randn();

    vn1 = exp(-P.k_gps*P.Ts_gps)*vn + P.sigma_n_gps * randn();
    ve1 = exp(-P.k_gps*P.Ts_gps)*ve + P.sigma_e_gps * randn();
    vh1 = exp(-P.k_gps*P.Ts_gps)*vh + P.sigma_h_gps * randn();
    vn = vn1;
    ve = ve1;
    vh = vh1;
    
    % construct total output
    y = [...
        y_gps_n;...
        y_gps_e;...
        y_gps_h;...
        y_gps_Vg;...
        y_gps_course;...
        ];
    
end



