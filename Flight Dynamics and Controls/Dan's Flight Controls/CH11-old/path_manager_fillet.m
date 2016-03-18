% path_manager_fillet
%   - follow lines between waypoints.  Smooth transition with fillets
%
%
% input is:
%   num_waypoints - number of waypoint configurations
%   waypoints    - an array of dimension 5 by P.size_waypoint_array.
%                - the first num_waypoints rows define waypoint
%                  configurations
%                - format for each waypoint configuration:
%                  [wn, we, wd, dont_care, Va_d]
%                  where the (wn, we, wd) is the NED position of the
%                  waypoint, and Va_d is the desired airspeed along the
%                  path.
%
% output is:
%   flag - if flag==1, follow waypoint path
%          if flag==2, follow orbit
%   
%   Va^d - desired airspeed
%   r    - inertial position of start of waypoint path
%   q    - unit vector that defines inertial direction of waypoint path
%   c    - center of orbit
%   rho  - radius of orbit
%   lambda = direction of orbit (+1 for CW, -1 for CCW)
%
function out = path_manager_fillet(in,P,start_of_simulation)

  NN = 0;
  num_waypoints = in(1+NN);
  waypoints = reshape(in(2+NN:5*P.size_waypoint_array+1+NN),5,P.size_waypoint_array);
  NN = NN + 1 + 5*P.size_waypoint_array;
  pn        = in(1+NN);
  pe        = in(2+NN);
  h         = in(3+NN);
  % Va      = in(4+NN);
  % alpha   = in(5+NN);
  % beta    = in(6+NN);
  % phi     = in(7+NN);
  % theta   = in(8+NN);
  % chi     = in(9+NN);
  % p       = in(10+NN);
  % q       = in(11+NN);
  % r       = in(12+NN);
  % Vg      = in(13+NN);
  % wn      = in(14+NN);
  % we      = in(15+NN);
  % psi     = in(16+NN);
  state     =  in(1+NN:16+NN);
  NN = NN + 16;
  t         = in(1+NN);
 
  p = [pn; pe; -h];
  
  R = 30; %Fillet Radius


  persistent waypoints_old   % stored copy of old waypoints
  persistent ptr_a           % waypoint pointer
  persistent state_transition % state of transition state machine
  persistent flag_need_new_waypoints % flag that request new waypoints from path planner
  
  
  if start_of_simulation || isempty(waypoints_old),
      waypoints_old = zeros(5,P.size_waypoint_array);
      flag_need_new_waypoints = 0;
     
  end
  
  % if the waypoints have changed, update the waypoint pointer
  if min(min(waypoints==waypoints_old))==0,
      ptr_a = 1;
      waypoints_old = waypoints;
      state_transition = 1;
      flag_need_new_waypoints = 0;
  end
  
  % define current and next two waypoints

  normalize = @(a, b) (a-b)/norm(a-b);
  q_i1 = normalize(waypoints(ptr_a), waypoints(ptr_a-1));
  q_i = normalize(waypoints(ptr_a+1), waypoints(ptr_a));
  script_phi = acos(-q_i1'*q_i);
  
  % define transition state machine
  switch state_transition,
      case 1, % follow straight line from wpp_a to wpp_b
          flag   = 1;  % following straight line path
          Va_d   = P.Va; % desired airspeed along waypoint path
          r      = waypoints(ptr_a-1);
          %q      = waypoints(ptr_a+1)-waypoints(ptr_a);
          %q      = q/norm(q);
          q = q_i1;
          c      = 0;
          rho    = inf;
          lambda = 0;
          z = waypoints(ptr_a)-(R/tan(script_phi/2))*q_i1;
          
             
      case 2, % follow orbit from wpp_a-wpp_b to wpp_b-wpp_c
          flag   = 2;  % following orbit
          Va_d   = P.Va; % desired airspeed along waypoint path
          r      = 0;
          %q      = ;
          %q      = q/norm(q);
          q = q_i1;
          q_next = q_i;
          %q_next = ;
          %q_next = q_next/norm(q_next);
          beta   = 0;
          c      = waypoints(ptr_a)+(R/sin(script_phi/2))*normalize(q_i1,q_i);
          rho    = R;
          lambda = sign(q_i1(1)*q_i(2)-q_i1(2)*q_i(1));
          z = waypoints(ptr_a)+(R/tan(script_phi/2))*q_i;
  end
  
  out = [flag; Va_d; r; q; c; rho; lambda; state; flag_need_new_waypoints];

end