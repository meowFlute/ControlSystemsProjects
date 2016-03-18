% sensors.m
%   Compute the output of rate gyros, accelerometers, and pressure sensors
%
%  Revised:
%   3/5/2010  - RB 
%   5/14/2010 - RB

function y = sensors(uu, P)

    % relabel the inputs
%    pn      = uu(1);
%    pe      = uu(2);
    pd      = uu(3);
%    u       = uu(4);
%    v       = uu(5);
%    w       = uu(6);
    phi     = uu(7);
    theta   = uu(8);
%    psi     = uu(9);
    p       = uu(10);
    q       = uu(11);
    r       = uu(12);
    F_x     = uu(13);
    F_y     = uu(14);
    F_z     = uu(15);
%    M_l     = uu(16);
%    M_m     = uu(17);
%    M_n     = uu(18);
    Va      = uu(19);
%    alpha   = uu(20);
%    beta    = uu(21);
%    wn      = uu(22);
%    we      = uu(23);
%    wd      = uu(24);
    
    % simulate rate gyros (units are rad/sec)
    y_gyro_x = sat(p + P.sigma_gyro*randn(), P.gyro_max, -P.gyro_max);
    y_gyro_y = sat(q + P.sigma_gyro*randn(), P.gyro_max, -P.gyro_max);
    y_gyro_z = sat(r + P.sigma_gyro*randn(), P.gyro_max, -P.gyro_max);

    % simulate accelerometers (units of g)
    y_accel_x = sat(F_x/P.mass + P.gravity*sin(theta)          + P.sigma_accel*randn(), P.accel_max, -P.accel_max);
    y_accel_y = sat(F_y/P.mass + P.gravity*sin(theta)*cos(phi) + P.sigma_accel*randn(), P.accel_max, -P.accel_max);
    y_accel_z = sat(F_z/P.mass + P.gravity*cos(theta)*cos(phi) + P.sigma_accel*randn(), P.accel_max, -P.accel_max);

    % simulate pressure sensors
    y_static_pres = 101.325 - (P.rho*P.gravity*(-pd))/1000 + P.beta_static + P.sigma_static*randn(); % kPa
    y_diff_pres   = (P.rho*(Va^2)/2)/1000 + P.beta_diff + P.sigma_diff*randn();                      % kPa

    % construct total output
    y = [...
        y_gyro_x;...
        y_gyro_y;...
        y_gyro_z;...
        y_accel_x;...
        y_accel_y;...
        y_accel_z;...
        y_static_pres;...
        y_diff_pres;...
    ];

end

function sat_value = sat(value,max,min)
    if value > max,
        value = max;
    end
    if value < min,
        value = min;
    end
    sat_value = value;
end

