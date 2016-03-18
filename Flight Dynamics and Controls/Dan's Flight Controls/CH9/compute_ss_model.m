function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(filename,x_trim,u_trim)
% x_trim is the trimmed state,
% u_trim is the trimmed input
  
% add stuff here  
[A,B,C,D] = linmod(filename,x_trim,u_trim);
%        pn pe  pd, u,  v,  w,  ph, th, ps, p,  q,  r
I_lat = [0, 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0;
         0, 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0;
         0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1;
         0, 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0;
         0, 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0];
I_lon = [0, 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0;
         0, 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0;
         0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0;
         0, 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0;
         0, 0,  -1,  0,  0,  0,  0,  0,  0,  0,  0,  0];
     
A_lon = I_lon*A*I_lon';
B_lon = I_lon*B;
A_lat = I_lat*A*I_lat';
B_lat = I_lat*B;