% path follow
%  - follow straight line path or orbit
%
% Modified:
%   3/25/2010  - RB
%   6/5/2010   - RB
%   11/08/2010 - RB
%   14/11/2014 - RWB
%
% input is:
%   flag - if flag==1, follow waypoint path
%          if flag==2, follow orbit
%   
%   Va^d   - desired airspeed
%   r      - inertial position of start of waypoint path
%   q      - unit vector that defines inertial direction of waypoint path
%   c      - center of orbit
%   rho    - radius of orbit
%   lambda - direction of orbit (+1 for CW, -1 for CCW)
%   xhat   - estimated MAV states (pn, pe, pd, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi)
%
% output is:
%  Va_c - airspeed command
%  h_c  - altitude command
%  chi_c - heading command
%  phi_ff - feed forward roll command
%
function out = path_follow(in,P)
  
  NN = 0;
  flag      = in(1+NN);
  Va_d      = in(2+NN);
  r_path    = [in(3+NN); in(4+NN); in(5+NN)];
  q_path    = [in(6+NN); in(7+NN); in(8+NN)];
  c_orbit   = [in(9+NN); in(10+NN); in(11+NN)];
  rho_orbit = in(12+NN);
  lam_orbit = in(13+NN);
  NN = NN + 13;
  pn        = in(1+NN);
  pe        = in(2+NN);
  h         = in(3+NN);
  Va        = in(4+NN);
  % alpha   = in(5+NN);
  % beta    = in(6+NN);
  phi       = in(7+NN);
  theta     = in(8+NN);
  chi       = in(9+NN);
  % p       = in(10+NN);
  % q       = in(11+NN);
   r       = in(12+NN);
  % Vg      = in(13+NN);
  % wn      = in(14+NN);
  % we      = in(15+NN);
  % psi     = in(16+NN);
  NN = NN + 16;
  t         = in(1+NN);
  
  switch flag,
      case 1, % follow straight line path specified by r and q
          %Tune these parameters
          chi_inf = 45 * pi/180;
          k_path = .1;
          
          
          p_path = [pn; pe; h];
          qn = q_path(1);
          qe = q_path(2);
          
          m = -2;
          chi_q = 0;
          while(((chi_q-chi)< -pi||(chi_q-chi)>pi)||chi_q ==0)
                m = m + 1;
                chi_q = atan2(qe,qn)+2*pi*m;
          end
          
          Rip = [cos(chi_q),    sin(chi_q), 0;...
                 -sin(chi_q),   cos(chi_q), 0;...
                 0,             0,          1];
          ep = Rip*(p_path-r_path);
          epy = ep(2);
          %epydot = Vg*sin(chi-chi_q);
          
          chi_c = chi_q-chi_inf*2/pi*atan(k_path*epy);
          h_c = r_path(3);
          phi_ff = phi;
           
      case 2, % follow orbit specified by c, rho, lam

          %Tune these parameters
          k_orbit = 3;
          
          lambda = lam_orbit;
          rho = rho_orbit;
          d = sqrt((pn-c_orbit(1))^2+(pe-c_orbit(2))^2);
          
          phaseAngle = atan2(pe-c_orbit(2),pn-c_orbit(1));
          while (phaseAngle-chi < pi)
             phaseAngle = phaseAngle + 2*pi; 
          end
          while (phaseAngle-chi > pi)
             phaseAngle = phaseAngle - 2*pi; 
          end
          chi_c = phaseAngle + lambda*(pi/2+atan(k_orbit*(d-rho)/rho));
          h_c = c_orbit(3); %h_c = c_d
          phi_ff = phi;
  end
  
  % command airspeed equal to desired airspeed
  Va_c = Va_d;
  
  % create output
  out = [Va_c; h_c; chi_c; phi_ff];
end


