function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(filename,x_trim,u_trim)
% x_trim is the trimmed state,
% u_trim is the trimmed input
  
[A,B,C,D] = linmod(filename, x_trim, u_trim);

% lateral linmod example from the book

% v
% p
% r
% phi
% psi

% delta_a
% delta_r

E1 = [...
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;...
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;];
E2 = [...
    0, 1, 0, 0;...
    0, 0, 1, 0;];
A_lat = E1 * A * E1';
B_lat = E1 * B * E2';

% longitudinal state space model 

% u
% w
% q
% theta
% h = -pd

% delta_e
% delta_t

E3 = [...
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;...
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
E4 = [...
    1, 0, 0, 0;...
    0, 0, 0, 1;];
A_lon = E3 * A * E3';
B_lon = E3 * B * E4';