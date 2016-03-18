function tau=arm_ctrl(in,P)
    theta_c = in(1);
    theta   = in(2);
    t       = in(3);
    
    
    %------------------------------------------------------------------
    % derivative control
    persistent thetadot
    persistent theta_d1
    % reset persistent variables at start of simulation
    if t<P.Ts,
        thetadot    = 0;
        theta_d1    = 0;
    end
    % update derivative of y
    thetadot = (2*P.tau-P.Ts)/(2*P.tau+P.Ts)*thetadot + 2/(2*P.tau+P.Ts)*(theta-theta_d1);
    theta_d1 = theta;

    
    %------------------------------------------------------------------
    % controller
    %
    % initialize registers to equilibrium values
    persistent cnum
    persistent cden
    if t<P.Ts,
        % compute equilibrium torque tau_e
        cnum = 0*ones(length(P.Cd_num),1);
        cden = 0*ones(length(P.Cd_den)-1,1);
    end
    % compute the error (note the use of the filtered command)
    error = theta_c-theta;
    % shift the control input register
    cnum = [error; cnum(1:end-1)];
    %
    % compute the output of the controller
    tau_e = P.m*P.g*(P.ell/2)*cos(theta_c);
    tau_tilde = P.Cd_num*cnum - P.Cd_den(2:end)*cden;
    tau =  -P.kd*thetadot + tau_e + tau_tilde;
    % anti-windup scheme:
    % reduce the error value stored in cnum(1), by exactly the
    % value needed to keep the torque at the limit.
    if tau_tilde>= P.tau_limit, 
        del = (P.tau_limit-tau_tilde)/(P.Cd_num(1)*cnum(1));
        cnum(1) = (1+del)*cnum(1);
        tau_tilde = P.tau_limit;
    end
    if tau_tilde<=-P.tau_limit,
        del = (-P.tau_limit-tau_tilde)/(P.Cd_num(1)*cnum(1));
        cnum(1) = (1+del)*cnum(1);
        tau_tilde = -P.tau_limit;
    end
    % shift the denominator register    
    cden = [tau_tilde; cden(1:end-1)];
    
end
