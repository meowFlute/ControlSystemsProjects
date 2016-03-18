function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(filename,x_trim,u_trim)
% x_trim is the trimmed state,
% u_trim is the trimmed input
  
% add stuff here
%x_trim(4) = -x_trim(4);
[A,B,C,D] = linmod(filename,x_trim,u_trim);

F = zeros(12,12);
count = 1;
for i=1:12
    for j=1:12
        F(i,j) = count;
        count = count + 1;
    end
end

%
I_lon = [0,  0,  0,  1, 0, 0, 0,   0,     0,   0, 0, 0;
         0,  0,  0,  0, 0, 1, 0,   0,     0,   0, 0, 0;
         0,  0,  0,  0, 0, 0, 0,   0,     0,   0, 1, 0;
         0,  0,  0,  0, 0, 0, 0,   1,     0,   0, 0, 0;
         0,  0,  1,  0, 0, 0, 0,   0,     0,   0, 0, 0];
%        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r
I_lat = [0,  0,  0,  0, 1, 0, 0,   0,     0,   0, 0, 0;
         0,  0,  0,  0, 0, 0, 0,   0,     0,   1, 0, 0;
         0,  0,  0,  0, 0, 0, 0,   0,     0,   0, 0, 1;
         0,  0,  0,  0, 0, 0, 1,   0,     0,   0, 0, 0;
         0,  0,  0,  0, 0, 0, 0,   0,     1,   0, 0, 0];
% pre-multiply pulls out rows
% post-multiply pulls out columns
A_lon = I_lon * A * I_lon';
B_lon = I_lon * B;
%B_lon = [B_lon(:,1),B_lon(:,4)];
    
A_lat = I_lat * A * I_lat';
B_lat = I_lat * B;
%}

end