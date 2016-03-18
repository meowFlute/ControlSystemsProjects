function out=arm_ctrl(in,P)
    theta_c  = in(1);
    theta    = in(2);
    t        = in(3);
    x        = in(4:5);
    
    xhat = x;
    
    tau_e  = P.m*P.g*P.ell/2*cos(theta_c);
    
    % add an observer
    persistent xhat_
    persistent P_
    persistent tau
    N = 10;
    y = theta;
    if t<P.Ts,
        xhat_  = [0; 0];
        P_ = .001*[1, 0; 0 1];
        tau = 0;
    else
        for i=1:N,
            xhat = xhat + P.Ts/N*f(xhat,tau,P);
            A = [0, 1; 3/2*P.g/P.ell*sin(xhat_(1)), -3*P.b/P.m/P.ell^2];
            P_ = P_ + P.Ts/N*(A*P_ + P_*A' + P.Q);
        end
        L = P_*P.C'/(P.C*P_*P.C'+P.R)
        P_ = (eye(2)-L*P.C)*P_
        xhat_ = xhat_ + L*(y-P.C*xhat_);
    end
   
    xhat_
    
    
    % define and initialize persistent variables
    persistent integrator
    persistent error_d1
    if t<P.Ts,
        integrator = 0;
        error_d1   = 0;
    end

    % compute the error and update the integrator
    error = theta_c-xhat(1);
    if abs(x(2))<.1,
        integrator = integrator + (P.Ts/2)*(error+error_d1);
    end
    error_d1 = error;
    
    % compute state feedback control with integrator
     tau_unsat = tau_e + P.kr*theta_c - P.K*xhat + P.ki*integrator;
    tau = sat(tau_unsat,P.tau_max);
    
    % integrator anti-windup
    if P.ki~=0,
        integrator = integrator + P.Ts/P.ki*(tau-tau_unsat);
    end
    
    out = [tau; xhat_];
    
end
    
function out = sat(in,limit)
    if     in > limit,      out = limit;
    elseif in < -limit,     out = -limit;
    else                    out = in;
    end
end

%----------------------------------------------------------
function xdot = f(x,u,P)

  theta    = x(1);
  thetadot = x(2);
  tau      = u(1);
  
  thetaddot = (3/P.m/P.ell^2)*(tau - P.b*thetadot - P.m*P.g*P.ell/2*cos(theta)); 

  xdot = [thetadot; thetaddot];

end

    
