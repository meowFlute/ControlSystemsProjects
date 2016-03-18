function [x_trim,u_trim,y_trim] = compute_trim(filename, Va, gamma, R)
% Va is the desired airspeed (m/s)
% gamma is the desired flight path angle (radians)
% R is the desired radius (m) - use (+) for right handed orbit, 
%                                   (-) for left handed orbit


% add stuff here

% Initial starting guess for states
x0 = [0;...   %Pn
      0;...   %Pe
      0;...   %Pd
      Va;...   %u  Approx: Airspeed should all be in U direction
      0;...    %v
      0;...    %w
      0;...    %phi
      gamma;...    %theta Approx: Assumes no wind and no angle of attack
      0;...    %psi
      0;...    %p
      0;...    %q
      0 ...     %r
      ];
  ix = [];  % Don't try to hold any states constant
  
  u0 = [0; 0; 0; 1];
  iu = []; %U values to hold constant
  
  y0 =[Va;   %airspeed
      gamma;    %alpha
      0];        %beta
  iy = [1;3];
  
  dx0 = [0;                 %Pn
         0;                 %Pe
         -Va*sin(gamma);    %Pd
         0;                 %u
         0;                 %v
         0;                 %w
         0;                 %phi
         0;                 %theta
         Va/R*cos(gamma);   %psi
         0;                 %p
         0;                 %q
         0                  %r
         ];
  idx = [3;4;5;6;7;8;9;10;11;12];

% compute trim conditions
[x_trim,u_trim,y_trim,dx_trim] = trim(filename,x0,u0,y0,ix,iu,iy,dx0,idx);

% check to make sure that the linearization worked (should be small)
norm(dx_trim(3:end)-dx0(3:end))

