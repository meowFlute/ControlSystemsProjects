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

    mass = P.mass;
    g = P.gravity;
    rho = P.rho;
    
    sigma_gyro = P.sigma_gyro
    % simulate rate gyros (units are rad/sec)
    y_gyro_x = p + sigma_gyro*randn();
    y_gyro_y = q + sigma_gyro*randn();
    y_gyro_z = r + sigma_gyro*randn();

    sigma_accel = P.sigma_accel;
    % simulate accelerometers (units of g)
    y_accel_x = F_x/mass+g*sin(theta)+sigma_accel*randn();
    y_accel_y = F_y/mass-g*cos(theta)*sin(phi)+sigma_accel*randn();
    y_accel_z = F_z/mass-g*cos(theta)*cos(phi)+sigma_accel*randn();

    beta_diff = .020;
    sigma_diff = P.sigma_diff;
    beta_abs = .125;
    sigma_abs = P.sigma_abs;
    % simulate pressure sensors
    y_static_pres = -rho*g*pd + beta_abs + sigma_abs*randn();
    y_diff_pres = rho*Va^2/2 + beta_diff + sigma_diff*randn();

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



