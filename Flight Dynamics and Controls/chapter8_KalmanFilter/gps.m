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
    
    persistent vn1;
    persistent ve1;
    persistent vh1;
    persistent time_dl;
 
    if(isempty(vn1))
       vn1 = 0;
       ve1 = 0;
       vh1 = 0;
    end
        
    kGPS = 1/1100;
    vn = exp(-kGPS*P.Ts_gps)*vn1 + 0.21*randn();
    ve = exp(-kGPS*P.Ts_gps)*ve1 + 0.21*randn();
    vh = exp(-kGPS*P.Ts_gps)*vh1 + 0.40*randn();
    
    vn1 = vn;
    ve1 = ve;
    vh1 = vh;
    
    % construct North, East, and altitude GPS measurements
    y_gps_n =  pn  + vn;
    y_gps_e =  pe  + ve; 
    y_gps_h = -pd  + vh; 
    
    % construct groundspeed and course measurements
    sigma_v = 0.0025;
    sigma_chi = 0.0025;
    y_gps_Vg     = sqrt((Va*cos(psi) + wn)^2 + (Va*sin(psi) + we)^2) + sigma_v*randn();
    y_gps_course = (atan2(Va*sin(psi)+we,Va*cos(psi)+wn)+sigma_chi*randn());

    % construct total output
    y = [...
        y_gps_n;...
        y_gps_e;...
        y_gps_h;...
        y_gps_Vg;...
        y_gps_course;...
        ];
    
end