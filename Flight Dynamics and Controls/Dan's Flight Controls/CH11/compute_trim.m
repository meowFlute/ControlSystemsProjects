function [x_trim,u_trim,y_trim] = compute_trim(filename, Va, gamma, R)
% Va is the desired airspeed (m/s)
% gamma is the desired flight path angle (radians)
% R is the desired radius (m) - use (+) for right handed orbit, 
%                                   (-) for left handed orbit

x0 = [0;        % pn
      0;        % pe
      0;        % h
      Va;       % u
      0;        % v
      0;        % w
      0;        % phi
      gamma;    % theta
      0;        % psi
      0;        % p
      0;        % q
      0];       % r
% values to hold constant
ix = [];

u0 = [0;        % delta_e
      0;        % delta_a
      0;        % delta_r
      1];       % delta_t
% values to hold constant
iu = [];
  
y0 = [Va; % airspeed
      gamma;  % alpha
      0]; % beta
% values to hold constant
iy = [1;3];

dx0 = [0;
       0;
       -Va*sin(gamma);
       0;
       0;
       0;
       0;
       0;
       Va/R*cos(gamma);
       0;
       0;
       0];
% values to hold constant
idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12];

% compute trim conditions
[x_trim,u_trim,y_trim,dx_trim] = trim(filename,x0,u0,y0,ix,iu,iy,dx0,idx);

% check to make sure that the linearization worked (should be small)
norm(dx_trim(3:end)-dx0(3:end))

