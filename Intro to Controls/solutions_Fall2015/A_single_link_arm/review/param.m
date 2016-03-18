%--------------------------------------------------------
% actual parameters - used in physics model in simulation
% actual system parameters
AP.m = 0.5;  % kg
AP.ell = 0.3; % m
AP.b = 0.01; % N m s
AP.g = 9.8; % m/s^2
% actual initial conditions
AP.theta0 = 0;
AP.thetadot0 = 0;
% input constraint
AP.tau_max = 1;

%--------------------------------------------------------
% control parameters - these are the parameters known to the controller,
% which may not be the actual parameters
P.m = AP.m*(0.9);  % kg
P.ell = AP.ell*(0.95); % m
P.b = AP.b*(1.1); % N m s
P.g = AP.g; % m/s^2
% initial conditions
P.theta0 = AP.theta0;
P.thetadot0 = AP.thetadot0;
% input constraint
P.tau_max = AP.tau_max;

% equalibrium torque
P.theta_e = 0;
P.tau_e = P.m*P.g*P.ell/3*cos(P.theta_e);

% sample rate of controller
P.Ts = 0.01;
% dirty derivative gain
P.tau = 0.05;

% PID design
A_th = 50*pi/180;
zeta = 0.707;
P.kp = (P.tau_max-P.tau_e)/A_th;
wn = sqrt(3*P.kp/P.m/P.ell^2);
P.kd = P.m*P.ell^2/3*2*zeta*wn-P.b;
P.ki = 0.1;

switch 3,
    case 1, %PID

    
    case 2, % loopshaping
        
    % loopshaping parameters
    % transfer function
    G = tf(3/P.m/P.ell^2,[1,3*P.b/P.m/P.ell^2,0]);
    % derivative controlled plant
    G_kd = tf(3/P.m/P.ell^2,[1,3*(P.b+P.kd)/P.m/P.ell^2,0]);
    C_pi = tf([P.kp,P.ki],[1,0]);
    GC_pid = series(C_pi,G_kd);

    P.tau_limit=(P.tau_max-P.tau_e)/A_th;

    %C = 1;
    C_total = series(C,C_pi);

    Cd = c2d(C_total,P.Ts,'tustin');
    [P.Cd_num,P.Cd_den]=tfdata(Cd,'v');

    % draw bode plots
    figure(4), clf
    bode(G)
    hold on
    grid on
    bode(GC_pid)
    bode(series(C,GC_pid))
    legend('G(s)','C_{pid}(s)G(s)','C(s)C_{pid}(s)G(s)')

        
    case 3, % observer-based control
        
        % state space model
        P.A = [0 1; 0 -3*P.b/P.m/(P.ell)^2];
        P.B = [0;3/P.m/(P.ell^2)];
        P.C = [1 0];
        P.D = 0;
        
        % check observability and controllability
        if rank(ctrb(P.A,P.B))~=2,
            disp('The system is not controllable')
        end
        if rank(obsv(P.A,P.C))~=2,
            disp('The system is not observable')
        end

        % augmented state space model to include integrator
        A_int = [P.A, zeros(2,1); P.C, 0];
        B_int = [P.B; 0];

        % desired closed loop pole locations
        wn = 6.5;
        zeta = 0.707;
        p_cl = roots([1,2*zeta*wn,wn^2]);
        p_int = -.1; % position of integrator pole
        tmp = place(A_int,B_int,[p_cl;p_int]);
        P.K = tmp(:,1:2);
        P.kr = -1/(P.C*inv(P.A-P.B*P.K)*P.B);
        P.ki = tmp(1,3);

        % observer
        p_ob = 10*p_cl;
        P.L = place(P.A',P.C',p_ob)';
        poles_of_controller=eig(P.A-P.B*P.K-P.L*P.C)

        % observer with disturbance estimate
        P.Ao = [P.A, P.B; zeros(1,3)];
        P.Bo = [P.B; 0];
        P.Co = [P.C, 0];
        P.Lo = place(P.Ao',P.Co',[p_ob;-10])';
end




