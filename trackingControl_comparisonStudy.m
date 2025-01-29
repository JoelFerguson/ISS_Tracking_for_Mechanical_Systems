%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ISS tracking control for mechanical systems
% Description: Comparison study code for the paper 'Input-to-State Stable 
% tracking control design for fully-actuated mechanical systems using 
% position measurements only',
% submitted to XXX
% Authours: Joel Ferguson, Naoki Sakata
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise
% Clear workspace
clear
close all
clc

% Set Figure default values
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultLineLineWidth',2.0);
set(0,'DefaultAxesLineWidth',0.5);
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaultAxesNextPlot','add')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Ferguson 2025 Tuning 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation settings
sim.isSave = 0;        % 0 = Dont save results; 1 = Save results
sim.disturbance = 0;   % 0 = No disturbance; 1 = Disturbance active

% Simulation initial conditions
sim.q0 = [0 0].';   % configuration (q)
sim.p0 = [-1 2].';  % momentum (p)
sim.ph0 = [0 0].';  % momentum estimate (\hat p)
sim.t_start = 0.01;
sim.t_end = 10;      % simulation length

% Observer parameters
obs.kappa = 1;

% Controller parameters
ctrl.alpha = 0.1;
ctrl.Kd = 10*eye(2);
ctrl.Kp = 6*eye(2);

% Target trajectory
ctrl.qd = @(t) [sin(t); cos(2*t)];

% Disturbance torque input
if sim.disturbance
    sys.delta_p0 = @(t) [1;1];
else
    sys.delta_p0 = @(t) [0; 0];
end

%% Define 2 degree-of-freedom manipulator model
% Define degrees of freedom
sys.N = length(sim.q0);

% symbolic variables
syms q_sym [sys.N,1]
syms p_sym [sys.N,1]

% Parameters
% kinematic parameters
sys.l1 = 1;
sys.l2 = 1;
sys.G0 = @(q) [1 -1; 0 1];

% inertial parameters
sys.m1 = 3;
sys.m2 = 3;
sys.j1 = 3/12;
sys.j2 = 3/12;
sys.M = @(q) [sys.j1 + (sys.l1^2*sys.m1)/4 + sys.l1^2*sys.m2,   (sys.l1*sys.l2*sys.m2*cos(q(1) - q(2)))/2;
                (sys.l1*sys.l2*sys.m2*cos(q(1) - q(2)))/2,                      (sys.m2*sys.l2^2)/4 + sys.j2];
            
% friction parameters
sys.dj1 = @(q) 1;
sys.dj2 = @(q) 1;
sys.D0 = @(q) [sys.dj1(q)+sys.dj2(q), -sys.dj2(q); -sys.dj2(q), sys.dj2(q)];

% Potential parameters
sys.g = 9.8;
sys.V = @(q) sys.g*sys.m2*(sys.l1*sin(q(1)) + (sys.l2*sin(q(2)))/2) + (sys.g*sys.l1*sys.m1*sin(q(1)))/2;

% System energy
sys.H = @(q,p) 0.5*p.'*(sys.M(q)\p) + sys.V(q);
sys.T = @(q,p) 0.5*p.'*(sys.M(q)\p);

% energy gradients
sys.dVdq = matlabFunction(jacobian(sys.V(q_sym),q_sym).','vars',{q_sym});
sys.dHdq = matlabFunction(jacobian(sys.H(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);
sys.dHdp = @(q,p) sys.M(q)\p;

% system ODE
sys.dx = @(t,q,p,u) [zeros(sys.N) eye(sys.N); -eye(sys.N) -sys.D0(q)]*[sys.dHdq(q,p); sys.dHdp(q,p)] + [zeros(sys.N); sys.G0(q)]*u + [zeros(sys.N,1); sys.delta_p0(t)];

% Define centrifugal/Coriolis matrix needed for EL methods
sys.C = @(q,dq) [0, 0.5*sys.l1*sys.l2*sys.m2*sin(q(1)-q(2))*dq(2);
                -0.5*sys.l1*sys.l2*sys.m2*sin(q(1)-q(2))*dq(1), 0];

%% Define momentum observer system
% Compute mass matrix square root
obs.Ti = @(q) sqrtm(sys.M(q));

% Compute mass matrix derivatives
obs.dMdqi = cell(sys.N,1);
obs.dTidqi = cell(sys.N,1);
for i=1:sys.N
    obs.dMdqi{i} = matlabFunction(diff(sys.M(q_sym),q_sym(i)),'vars',{q_sym});
    obs.dTidqi{i} = @(q) lyap(obs.Ti(q), -feval(obs.dMdqi{i},q));
end

% Function handle for A matrix
obs.A = @(q,p) constructAmatrix(q,p,sys.N,obs.dTidqi);

% Compute S matrix
obs.S = @(q,p) (obs.Ti(q)\(obs.A(q,p).' - obs.A(q,p)))/obs.Ti(q);

% Function handle for \bar S matrix
obs.Sb = @(q,ph) constructSbMatrix(q,ph,sys.N,obs.S);

% Compute remaining system matrices in transformed coordinates
obs.D = @(q) (obs.Ti(q)\sys.D0(q))/obs.Ti(q);
obs.G = @(q) obs.Ti(q)\sys.G0(q);

% Observer energy function
obs.ph = @(q,xp,phi) xp + phi*q;
obs.Hp = @(p,ph) 0.5*(p-ph).'*(p-ph);

% Condition for switching
obs.switchCond = @(q,ph,phi) min(eig(phi*eye(2)/obs.Ti(q) - 0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).'))) - obs.kappa;
obs.stopEventWrapper = @(t,x) stopEvent(t,x,obs.switchCond,obs.ph);

% Observer dynamics
obs.dxp = @(q,ph,u,uo,phi) (obs.S(q,ph) - obs.D(q) - phi*eye(2)/obs.Ti(q))*ph - obs.Ti(q)\sys.dVdq(q) + obs.G(q)*(u + uo);
obs.dphi = 0;

% compute initial phi
sim.phi0 = 0;
while obs.switchCond(sim.q0,sim.ph0,sim.phi0) <= 0
    sim.phi0 = sim.phi0 + obs.kappa;
end

% Define observer input
obs.uo = @(q,ye) -obs.G(q)\ye;

%% Define the tracking controller
% symbolic variables
syms t_sym

% Compute the trajectory derivatives
ctrl.dqd = matlabFunction(jacobian(ctrl.qd(t_sym),t_sym),'vars',t_sym);                 % Derivative of qd with respect to t
ctrl.ddqd = matlabFunction(jacobian(ctrl.dqd(t_sym),t_sym),'vars',t_sym);               % Second derivative of qd with respect to t

% Compute the desired momentum as a function of the desired velocity
ctrl.pd = @(q,t) obs.Ti(q)*ctrl.dqd(t);                                                 % Desired momentum (pd)
ctrl.dpddq = matlabFunction(jacobian(ctrl.pd(q_sym,t_sym),q_sym),'vars',{q_sym,t_sym}); % Derivative of pd with respect to q

% Define tracking error coordinates
ctrl.qe = @(q,t) q - ctrl.qd(t);
ctrl.pe = @(q,ph,t) ph - ctrl.pd(q,t);

% Define feed-forward tracking control law
ctrl.u = @(uo,v,q,ph,t) -uo + obs.G(q)\(-(obs.S(q,ph)-obs.D(q))*ctrl.pd(q,t) + obs.Ti(q)\sys.dVdq(q) + ctrl.dpddq(q,t)*(obs.Ti(q)\ph) + obs.Ti(q)*ctrl.ddqd(t) + v);

% Define KPES control law
ctrl.v = @(q,ph,t) ctrl.alpha*(obs.S(q,ph)-obs.D(q)-ctrl.Kd)*ctrl.Kp*(ctrl.qe(q,t) + ctrl.alpha*ctrl.pe(q,ph,t)) - obs.Ti(q)\ctrl.Kp*(ctrl.qe(q,t) + ctrl.alpha*ctrl.pe(q,ph,t)) - ctrl.Kd*ctrl.pe(q,ph,t);

% Define tracking control storage function
ctrl.He = @(qe, pe) 0.5*pe.'*pe + 0.5*(qe + ctrl.alpha*pe).'*ctrl.Kp*(qe + ctrl.alpha*pe);
ctrl.dHedqe = @(qe,pe) ctrl.Kp*(qe + ctrl.alpha*pe);
ctrl.dHedpe = @(qe,pe) pe + ctrl.alpha*ctrl.Kp*(qe + ctrl.alpha*pe);

% Define passive output from tracking error dynamics
crtl.Q = [ctrl.Kp ctrl.alpha*ctrl.Kp;
        ctrl.alpha*ctrl.Kp eye(sys.N)+ctrl.alpha^2*ctrl.Kp];
ctrl.ye = @(q,qe,pe,ph,phi,t) -obs.Ti(q)\[eye(sys.N), phi*eye(sys.N) - ctrl.dpddq(q,t).']*[ctrl.dHedqe(qe,pe); ctrl.dHedpe(qe,pe)];

%% Run simulation
% initial conditions (see top of script)
sim.xp0 = sim.ph0 - sim.phi0*sim.q0;
sim.x0 = [sim.q0; sim.p0; sim.xp0; sim.phi0];

% Resolve inputs to be a function of states for ODE solver
sim.uo = @(q,xp,phi,t) obs.uo(q,ctrl.ye(q,ctrl.qe(q,t),ctrl.pe(q,obs.ph(q,xp,phi),t),obs.ph(q,xp,phi),phi,t));
sim.v = @(q,xp,phi,t) ctrl.v(q,obs.ph(q,xp,phi),t);
sim.u = @(q,xp,phi,t) ctrl.u(sim.uo(q,xp,phi,t),sim.v(q,xp,phi,t),q,obs.ph(q,xp,phi),t);

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),sim.u(x(1:2),x(5:6),x(7),t));
                obs.dxp(x(1:2),obs.ph(x(1:2),x(5:6),x(7)),sim.u(x(1:2),x(5:6),x(7),t),sim.uo(x(1:2),x(5:6),x(7),t),x(7));
                obs.dphi];

% stop solver when gain update is required
sim.options = odeset('RelTol',1e-6,'Events',obs.stopEventWrapper);
% Start solver
[t,x] = ode45(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);
% Store outputs
x_cat = x;
t_cat = t;
% Iteratively continue solve until final time point reached
while(t(end) < sim.t_end)
    % Extract initial conditions for next solve
    x_end = x(end,:).';
    q_end = x_end(1:sys.N);
    p_end = x_end(sys.N+1:2*sys.N);
    xp_end = x_end(2*sys.N+1:3*sys.N);
    phi_end = x_end(3*sys.N+1);
    ph_end = obs.ph(q_end,xp_end,phi_end);

    if norm(x_end)>100
        disp("State is deverging. Stop ode")
        return
    end
    
    % update observer values
    phi_minus = phi_end;
    phi_plus = phi_minus + obs.kappa;
    xp_plus = xp_end - obs.kappa*q_end;
    
    % Construct new initial state vector
    sim.x0 = [q_end; p_end; xp_plus; phi_plus];

    % stop solver when gain update is required
    sim.options = odeset('RelTol',1e-6,'Events',obs.stopEventWrapper);
    % Start solver
    [t,x] = ode45(sim.ode,[t(end) sim.t_end],sim.x0,sim.options);

    % Store outputs
    x_cat = [x_cat; x];
    t_cat = [t_cat; t];
end
t = t_cat;
x = x_cat;

%% Plot output
% Unpack state vector
resFerguson.t = t;
resFerguson.q = x(:,1:sys.N);
resFerguson.p0 = x(:,sys.N+1:2*sys.N);
resFerguson.xp = x(:,2*sys.N+1:3*sys.N);
resFerguson.phi = x(:,3*sys.N+1);

% Preallocate space for results
resFerguson.H = zeros(length(t),1);
resFerguson.Ho = zeros(length(t),1);
resFerguson.He = zeros(length(t),1);
resFerguson.W = zeros(length(t),1);
resFerguson.p0h = zeros(length(t),2);
resFerguson.p = zeros(length(t),2);
resFerguson.ph = zeros(length(t),2);
resFerguson.eig_Sb = zeros(length(t),1);
resFerguson.switchCond = zeros(length(t),1);
resFerguson.qe = zeros(length(t),2);
resFerguson.pe = zeros(length(t),2);
resFerguson.qd = zeros(length(t),2);
resFerguson.qe_norm = zeros(length(t),1);
resFerguson.uNorm = zeros(length(t),1);
resFerguson.uRef = zeros(length(t),2);

% Evaluate functions of interest for each time step
for i=1:length(t)
    % Unpack state for the ith time step
    t = resFerguson.t(i);
    q = resFerguson.q(i,:).';
    p0 = resFerguson.p0(i,:).';
    xp = resFerguson.xp(i,:).';
    phi = resFerguson.phi(i);

    % Compute coordinate transforms and error coordinates
    p = obs.Ti(q)\p0;
    ph = obs.ph(q,xp,phi);

    % Evaluate functions of interest
    resFerguson.H(i) = sys.H(q,p0);
    resFerguson.p(i,:) = p.';
    resFerguson.ph(i,:) = ph.';
    resFerguson.p0h(i,:) = (obs.Ti(q)*ph).';
    resFerguson.Ho(i) = obs.Hp(p,ph);
    resFerguson.u(i,:) = sim.u(q,xp,phi,t).';
    resFerguson.uNorm(i) = norm(resFerguson.u(i,:));
    resFerguson.qe(i,:) = ctrl.qe(q,t).';
    resFerguson.pe(i,:) = ctrl.pe(q,ph,t).';
    resFerguson.He(i) = ctrl.He(resFerguson.qe(i,:).', resFerguson.pe(i,:).');
    resFerguson.W(i) = resFerguson.Ho(i) + resFerguson.He(i);
    resFerguson.qd(i,:) = ctrl.qd(t).';
    resFerguson.qe_norm(i) = norm(resFerguson.qe(i,:));
    resFerguson.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resFerguson.uDiffNorm(i) = norm(resFerguson.u(i,:) - resFerguson.uRef(i,:));

    % compute max eigenvalues of symm[Sb]
    resFerguson.eig_Sb(i) = max(eig(0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).')));
    resFerguson.switchCond(i) = obs.switchCond(q,ph,phi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Ferguson 2025 Tuning 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define tuning gain nuber 2
% Observer parameters
obs.kappa = 0.3;

% Controller parameters
ctrl.alpha = 0.18;
ctrl.Kd = 0.3*eye(2);
ctrl.Kp = 1.4*eye(2);

%% Define momentum observer system
% Compute mass matrix square root
obs.Ti = @(q) sqrtm(sys.M(q));

% Compute mass matrix derivatives
obs.dMdqi = cell(sys.N,1);
obs.dTidqi = cell(sys.N,1);
for i=1:sys.N
    obs.dMdqi{i} = matlabFunction(diff(sys.M(q_sym),q_sym(i)),'vars',{q_sym});
    obs.dTidqi{i} = @(q) lyap(obs.Ti(q), -feval(obs.dMdqi{i},q));
end

% Function handle for A matrix
obs.A = @(q,p) constructAmatrix(q,p,sys.N,obs.dTidqi);

% Compute S matrix
obs.S = @(q,p) (obs.Ti(q)\(obs.A(q,p).' - obs.A(q,p)))/obs.Ti(q);

% Function handle for \bar S matrix
obs.Sb = @(q,ph) constructSbMatrix(q,ph,sys.N,obs.S);

% Compute remaining system matrices in transformed coordinates
obs.D = @(q) (obs.Ti(q)\sys.D0(q))/obs.Ti(q);
obs.G = @(q) obs.Ti(q)\sys.G0(q);

% Observer energy function
obs.ph = @(q,xp,phi) xp + phi*q;
obs.Hp = @(p,ph) 0.5*(p-ph).'*(p-ph);

% Condition for switching
obs.switchCond = @(q,ph,phi) min(eig(phi*eye(2)/obs.Ti(q) - 0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).'))) - obs.kappa;
obs.stopEventWrapper = @(t,x) stopEvent(t,x,obs.switchCond,obs.ph);

% Observer dynamics
obs.dxp = @(q,ph,u,uo,phi) (obs.S(q,ph) - obs.D(q) - phi*eye(2)/obs.Ti(q))*ph - obs.Ti(q)\sys.dVdq(q) + obs.G(q)*(u + uo);
obs.dphi = 0;

% compute initial phi
sim.phi0 = 0;
while obs.switchCond(sim.q0,sim.ph0,sim.phi0) <= 0
    sim.phi0 = sim.phi0 + obs.kappa;
end

% Define observer input
obs.uo = @(q,ye) -obs.G(q)\ye;

%% Define the tracking controller
% symbolic variables
syms t_sym

% Compute the trajectory derivatives
ctrl.dqd = matlabFunction(jacobian(ctrl.qd(t_sym),t_sym),'vars',t_sym);                 % Derivative of qd with respect to t
ctrl.ddqd = matlabFunction(jacobian(ctrl.dqd(t_sym),t_sym),'vars',t_sym);               % Second derivative of qd with respect to t

% Compute the desired momentum as a function of the desired velocity
ctrl.pd = @(q,t) obs.Ti(q)*ctrl.dqd(t);                                                 % Desired momentum (pd)
ctrl.dpddq = matlabFunction(jacobian(ctrl.pd(q_sym,t_sym),q_sym),'vars',{q_sym,t_sym}); % Derivative of pd with respect to q

% Define tracking error coordinates
ctrl.qe = @(q,t) q - ctrl.qd(t);
ctrl.pe = @(q,ph,t) ph - ctrl.pd(q,t);

% Define feed-forward tracking control law
ctrl.u = @(uo,v,q,ph,t) -uo + obs.G(q)\(-(obs.S(q,ph)-obs.D(q))*ctrl.pd(q,t) + obs.Ti(q)\sys.dVdq(q) + ctrl.dpddq(q,t)*(obs.Ti(q)\ph) + obs.Ti(q)*ctrl.ddqd(t) + v);

% Define KPES control law
ctrl.v = @(q,ph,t) ctrl.alpha*(obs.S(q,ph)-obs.D(q)-ctrl.Kd)*ctrl.Kp*(ctrl.qe(q,t) + ctrl.alpha*ctrl.pe(q,ph,t)) - obs.Ti(q)\ctrl.Kp*(ctrl.qe(q,t) + ctrl.alpha*ctrl.pe(q,ph,t)) - ctrl.Kd*ctrl.pe(q,ph,t);

% Define tracking control storage function
ctrl.He = @(qe, pe) 0.5*pe.'*pe + 0.5*(qe + ctrl.alpha*pe).'*ctrl.Kp*(qe + ctrl.alpha*pe);
ctrl.dHedqe = @(qe,pe) ctrl.Kp*(qe + ctrl.alpha*pe);
ctrl.dHedpe = @(qe,pe) pe + ctrl.alpha*ctrl.Kp*(qe + ctrl.alpha*pe);

% Define passive output from tracking error dynamics
crtl.Q = [ctrl.Kp ctrl.alpha*ctrl.Kp;
        ctrl.alpha*ctrl.Kp eye(sys.N)+ctrl.alpha^2*ctrl.Kp];
ctrl.ye = @(q,qe,pe,ph,phi,t) -obs.Ti(q)\[eye(sys.N), phi*eye(sys.N) - ctrl.dpddq(q,t).']*[ctrl.dHedqe(qe,pe); ctrl.dHedpe(qe,pe)];

%% Run simulation
% initial conditions (see top of script)
sim.xp0 = sim.ph0 - sim.phi0*sim.q0;
sim.x0 = [sim.q0; sim.p0; sim.xp0; sim.phi0];

% Resolve inputs to be a function of states for ODE solver
sim.uo = @(q,xp,phi,t) obs.uo(q,ctrl.ye(q,ctrl.qe(q,t),ctrl.pe(q,obs.ph(q,xp,phi),t),obs.ph(q,xp,phi),phi,t));
sim.v = @(q,xp,phi,t) ctrl.v(q,obs.ph(q,xp,phi),t);
sim.u = @(q,xp,phi,t) ctrl.u(sim.uo(q,xp,phi,t),sim.v(q,xp,phi,t),q,obs.ph(q,xp,phi),t);

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),sim.u(x(1:2),x(5:6),x(7),t));
                obs.dxp(x(1:2),obs.ph(x(1:2),x(5:6),x(7)),sim.u(x(1:2),x(5:6),x(7),t),sim.uo(x(1:2),x(5:6),x(7),t),x(7));
                obs.dphi];

% stop solver when gain update is required
sim.options = odeset('RelTol',1e-6,'Events',obs.stopEventWrapper);
% Start solver
[t,x] = ode45(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);
% Store outputs
x_cat = x;
t_cat = t;
% Iteratively continue solve until final time point reached
while(t(end) < sim.t_end)
    % Extract initial conditions for next solve
    x_end = x(end,:).';
    q_end = x_end(1:sys.N);
    p_end = x_end(sys.N+1:2*sys.N);
    xp_end = x_end(2*sys.N+1:3*sys.N);
    phi_end = x_end(3*sys.N+1);
    ph_end = obs.ph(q_end,xp_end,phi_end);

    if norm(x_end)>100
        disp("State is deverging. Stop ode")
        return
    end
    
    % update observer values
    phi_minus = phi_end;
    phi_plus = phi_minus + obs.kappa;
    xp_plus = xp_end - obs.kappa*q_end;
    
    % Construct new initial state vector
    sim.x0 = [q_end; p_end; xp_plus; phi_plus];

    % stop solver when gain update is required
    sim.options = odeset('RelTol',1e-6,'Events',obs.stopEventWrapper);
    % Start solver
    [t,x] = ode45(sim.ode,[t(end) sim.t_end],sim.x0,sim.options);

    % Store outputs
    x_cat = [x_cat; x];
    t_cat = [t_cat; t];
end
t = t_cat;
x = x_cat;

%% Plot output
% Unpack state vector
resFerguson2.t = t;
resFerguson2.q = x(:,1:sys.N);
resFerguson2.p0 = x(:,sys.N+1:2*sys.N);
resFerguson2.xp = x(:,2*sys.N+1:3*sys.N);
resFerguson2.phi = x(:,3*sys.N+1);

% Preallocate space for results
resFerguson2.H = zeros(length(t),1);
resFerguson2.Ho = zeros(length(t),1);
resFerguson2.He = zeros(length(t),1);
resFerguson2.W = zeros(length(t),1);
resFerguson2.p0h = zeros(length(t),2);
resFerguson2.p = zeros(length(t),2);
resFerguson2.ph = zeros(length(t),2);
resFerguson2.eig_Sb = zeros(length(t),1);
resFerguson2.switchCond = zeros(length(t),1);
resFerguson2.qe = zeros(length(t),2);
resFerguson2.pe = zeros(length(t),2);
resFerguson2.qd = zeros(length(t),2);
resFerguson2.qe_norm = zeros(length(t),1);
resFerguson2.uNorm = zeros(length(t),1);
resFerguson2.uRef = zeros(length(t),2);

% Evaluate functions of interest for each time step
for i=1:length(t)
    % Unpack state for the ith time step
    t = resFerguson2.t(i);
    q = resFerguson2.q(i,:).';
    p0 = resFerguson2.p0(i,:).';
    xp = resFerguson2.xp(i,:).';
    phi = resFerguson2.phi(i);

    % Compute coordinate transforms and error coordinates
    p = obs.Ti(q)\p0;
    ph = obs.ph(q,xp,phi);

    % Evaluate functions of interest
    resFerguson2.H(i) = sys.H(q,p0);
    resFerguson2.p(i,:) = p.';
    resFerguson2.ph(i,:) = ph.';
    resFerguson2.p0h(i,:) = (obs.Ti(q)*ph).';
    resFerguson2.Ho(i) = obs.Hp(p,ph);
    resFerguson2.u(i,:) = sim.u(q,xp,phi,t).';
    resFerguson2.uNorm(i) = norm(resFerguson2.u(i,:));
    resFerguson2.qe(i,:) = ctrl.qe(q,t).';
    resFerguson2.pe(i,:) = ctrl.pe(q,ph,t).';
    resFerguson2.He(i) = ctrl.He(resFerguson2.qe(i,:).', resFerguson2.pe(i,:).');
    resFerguson2.W(i) = resFerguson2.Ho(i) + resFerguson2.He(i);
    resFerguson2.qd(i,:) = ctrl.qd(t).';
    resFerguson2.qe_norm(i) = norm(resFerguson2.qe(i,:));
    resFerguson2.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resFerguson2.uDiffNorm(i) = norm(resFerguson2.u(i,:) - resFerguson2.uRef(i,:));

    % compute max eigenvalues of symm[Sb]
    resFerguson2.eig_Sb(i) = max(eig(0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).')));
    resFerguson2.switchCond(i) = obs.switchCond(q,ph,phi);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loria 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define controller parameters
ctrl.a = 1.5;
ctrl.b = 2;
ctrl.kd = 16;
ctrl.kc = 1;
ctrl.kp = 1;
ctrl.kdelta = 1;

% Define tracking error
ctrl.qt = @(t,q) q - ctrl.qd(t);

% Define dirty derivative
ctrl.dqc = @(t,q,qc) -ctrl.a*(qc + ctrl.b*ctrl.qt(t,q));
ctrl.dqc0 = [0;0];

% Define control input
ctrl.u = @(t,q,qc) sys.G0(q)\(-ctrl.kp*ctrl.qt(t,q) - ctrl.kd*(qc + ctrl.b*ctrl.qt(t,q)) + sys.M(q)*ctrl.ddqd(t) + sys.C(q,ctrl.dqd(t))*ctrl.dqd(t) + sys.dVdq(q));

%% Run simulation
% initial conditions (see top of script)
sim.x0 = [sim.q0; sim.p0; ctrl.dqc0];

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),ctrl.u(t,x(1:2),x(5:6)));
                ctrl.dqc(t,x(1:2),x(5:6))];

% Set ODE options
sim.options = odeset('RelTol',1e-6);
% Start solver
[t,x] = ode45(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);

% Unpack state vector
resLoria.t = t;
resLoria.q = x(:,1:sys.N);
resLoria.p0 = x(:,sys.N+1:2*sys.N);
resLoria.qc = x(:,2*sys.N+1:3*sys.N);

for i=1:length(resLoria.t)
    t = resLoria.t(i).';
    q = resLoria.q(i,:).';
    qc = resLoria.qc(i,:).';
    qd = ctrl.qd(t);

    resLoria.qd(i,:) = qd.';
    resLoria.qe_norm(i) = norm(q - qd);
    resLoria.u(i,:) = ctrl.u(t,q,qc);
    resLoria.uNorm(i) = norm(resLoria.u(i,:));
    resLoria.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resLoria.uDiffNorm(i) = norm(resLoria.u(i,:) - resLoria.uRef(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cruz 2021  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the tracking controller from Cruz 2021
% Define controller parameters
ctrl.r1 = 1.5;
ctrl.power = (2-ctrl.r1)/ctrl.r1;
% Set large bound so its not hit. Assumping the same bound for all inputs.
ctrl.delta = 100;
ctrl.P = diag([7, 4]);
ctrl.D = diag([4, 2]);
ctrl.B = eye(2);
ctrl.A = 7*eye(2);

% Define saturation function wrapper
sat_delta = @(x,p) sat_vec(x,p,ctrl.delta);

% Define tracking error
ctrl.qt = @(t,q) q - ctrl.qd(t);

% Filter dynamics
ctrl.nu = @(t,q,x) x + ctrl.B*ctrl.qt(t,q);
ctrl.dx = @(t,q,x) -ctrl.A*sat_delta(ctrl.nu(t,q,x),1/ctrl.r1);

% Define control input
ctrl.tau_ff = @(t,q) sys.dVdq(q) + sys.M(q)*ctrl.ddqd(t) + sys.C(q,ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(q)*ctrl.dqd(t);
ctrl.tau_pd = @(t,q,nu) -ctrl.P*sat_delta(ctrl.qt(t,q),ctrl.power) - ctrl.D*sat_delta(nu,ctrl.power);
ctrl.u = @(t,q,x) sys.G0(q)\(ctrl.tau_ff(t,q) + ctrl.tau_pd(t,q,ctrl.nu(t,q,x)));

% Control inital conditions
ctrl.x0 = [0;0];

%% Run simulation
% initial conditions (see top of script)
sim.x0 = [sim.q0; sim.p0; ctrl.x0];

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),ctrl.u(t,x(1:2),x(5:6)));
                ctrl.dx(t,x(1:2), x(5:6))];

% stop solver when gain update is required
sim.options = odeset('RelTol',1e-6);
% Start solver
[t,x] = ode45(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);

%% Plot output
% Unpack state vector
resCruz2021.t = t;
resCruz2021.q = x(:,1:sys.N);
resCruz2021.p0 = x(:,sys.N+1:2*sys.N);
resCruz2021.qc = x(:,2*sys.N+1:3*sys.N);

for i=1:length(resCruz2021.t)
    t = resCruz2021.t(i).';
    q = resCruz2021.q(i,:).';
    qc = resCruz2021.qc(i,:).';
    qd = ctrl.qd(t);

    resCruz2021.qd(i,:) = qd.';
    resCruz2021.qe_norm(i) = norm(q - qd);
    resCruz2021.u(i,:) = ctrl.u(t,q,qc);
    resCruz2021.uNorm(i) = norm(resCruz2021.u(i,:));
    resCruz2021.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resCruz2021.uDiffNorm(i) = norm(resCruz2021.u(i,:) - resCruz2021.uRef(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rascon 2020 IET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the tracking controller from Rascon2020IET
% Define controller parameters
ctrl.K1 = 50*eye(2);
ctrl.K2 = 10*eye(2);

% Define tracking error
ctrl.z1 = @(t,q) q - ctrl.qd(t);

% Define control input
ctrl.u = @(t,q,tau_c) sys.G0(q)\(sys.dVdq(q) + sys.M(q)*ctrl.ddqd(t) + sys.C(q,ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(q)*ctrl.dqd(t) + tau_c);
ctrl.tau_c = @(z1, z3) -ctrl.K1*z1 - z3;
crtl.dz3 = @(z1, z3) -ctrl.K1*z1 - ctrl.K2*z3;

% Total control law
ctrl.u_cat = @(t,q,z3) ctrl.u(t,q,ctrl.tau_c(ctrl.z1(t,q), z3));

% Control inital conditions
ctrl.z30 = [0;0];

%% Run simulation
% initial conditions (see top of script)
sim.x0 = [sim.q0; sim.p0; ctrl.z30];

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),ctrl.u_cat(t,x(1:2),x(5:6)));
                crtl.dz3(ctrl.z1(t,x(1:2)), x(5:6))];

% stop solver when gain update is required
sim.options = odeset('RelTol',1e-6);
% Start solver
[t,x] = ode45(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);

%% Plot output
% Unpack state vector
resRascon2020IET.t = t;
resRascon2020IET.q = x(:,1:sys.N);
resRascon2020IET.p0 = x(:,sys.N+1:2*sys.N);
resRascon2020IET.z3 = x(:,2*sys.N+1:3*sys.N);

for i=1:length(resRascon2020IET.t)
    t = resRascon2020IET.t(i).';
    q = resRascon2020IET.q(i,:).';
    qd = ctrl.qd(t);
    z1 = ctrl.z1(t,q);
    z3 = resRascon2020IET.z3(i,:).';
    tau_c = ctrl.tau_c(z1, z3);

    resRascon2020IET.qd(i,:) = qd.';
    resRascon2020IET.qe_norm(i) = norm(q - qd);
    resRascon2020IET.u(i,:) = ctrl.u(t,q,tau_c);
    resRascon2020IET.uNorm(i) = norm(resRascon2020IET.u(i,:));
    resRascon2020IET.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resRascon2020IET.uDiffNorm(i) = norm(resRascon2020IET.u(i,:) - resRascon2020IET.uRef(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rascon 2020 ISA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the tracking controller from Rascon 2020ISA
% Define controller parameters
ctrl.K1 = 10*eye(2);
ctrl.K2 = 10*eye(2);

% Define sign function
ctrl.nu = @(z1) [sign(z1(1)); sign(z1(2))];

% Define tracking error
ctrl.z1 = @(t,q) q - ctrl.qd(t);

% Define control input
ctrl.u = @(t,q,tau_c) sys.G0(q)\(sys.dVdq(q) + sys.M(q)*ctrl.ddqd(t) + tau_c);
ctrl.tau_c = @(t, q, z1) sys.C(q,ctrl.dqd(t))*ctrl.dqd(t) - ctrl.K1*ctrl.nu(z1) - ctrl.K2*z1;
crtl.dz3 = @(z1, z3) -ctrl.K1*z1;

% Total control law
ctrl.u_cat = @(t,q) ctrl.u(t,q,ctrl.tau_c(t, q, ctrl.z1(t,q)));

% Control inital conditions
ctrl.z30 = [0;0];

%% Run simulation
% initial conditions (see top of script)
sim.x0 = [sim.q0; sim.p0];

% Construct overall ode --- x = (q,p0,xp,phi)
sim.ode = @(t,x) [sys.dx(t,x(1:2),x(3:4),ctrl.u_cat(t,x(1:2)))];

% stop solver when gain update is required
sim.options = odeset('RelTol',1e-6); 
% Start solver
[res.t,res.x] = ode23s(sim.ode,[sim.t_start sim.t_end],sim.x0,sim.options);

%% Plot output
% Unpack state vector
resRascon2020ISA.t = res.t;
resRascon2020ISA.q = res.x(:,1:sys.N);
resRascon2020ISA.p0 = res.x(:,sys.N+1:2*sys.N);

for i=1:length(resRascon2020ISA.t)
    % Unpack state for the ith time step
    t = resRascon2020ISA.t(i);
    q = resRascon2020ISA.q(i,:).';
    resRascon2020ISA.u(i,:) = ctrl.u_cat(t,q).';
    resRascon2020ISA.qd(i,:) = ctrl.qd(t).';

    %resRascon2020ISA.qd(i,:) = qd.';
    resRascon2020ISA.qe_norm(i) = norm(resRascon2020ISA.q(i,:) - resRascon2020ISA.qd(i,:));
    resRascon2020ISA.u(i,:) = ctrl.u_cat(t,q).';
    resRascon2020ISA.uNorm(i) = norm(resRascon2020ISA.u(i,:));
    resRascon2020ISA.uRef(i,:) = (sys.G0(ctrl.qd(t))\(sys.M(ctrl.qd(t))*ctrl.ddqd(t) + sys.C(ctrl.qd(t),ctrl.dqd(t))*ctrl.dqd(t) + sys.D0(ctrl.qd(t))*ctrl.dqd(t) + sys.dVdq(ctrl.qd(t)))).';
    resRascon2020ISA.uDiffNorm(i) = norm(resRascon2020ISA.u(i,:) - resRascon2020ISA.uRef(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Plot comparison results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results
fig1 = figure("Name","Comparison plots");
fig1sub1 = subplot(1,2,1)
plot(resFerguson.t,resFerguson.uDiffNorm, resFerguson2.t,resFerguson2.uDiffNorm, resCruz2021.t,resCruz2021.uDiffNorm, resRascon2020IET.t,resRascon2020IET.uDiffNorm, resRascon2020ISA.t,resRascon2020ISA.uDiffNorm)
fig1sub1.YScale = 'log'
legend('Ferguson 2025 Tuning 1','Ferguson 2025 Tuning 2','Cruz 2021','Rascon 2020 IET','Rascon 2020 ISA','Interpreter','Latex','Orientation','vertical')
xlabel('time (s)')
ylabel('$$||u(t) - u_d(t)||$$','Interpreter','Latex')

% create smaller axes in top right, and plot on it
boxplot1 = axes('Position',[.17 .17 .1 .2])
box on
plot(resFerguson.t,resFerguson.uDiffNorm, resFerguson2.t,resFerguson2.uDiffNorm, resCruz2021.t,resCruz2021.uDiffNorm, resRascon2020IET.t,resRascon2020IET.uDiffNorm, resRascon2020ISA.t,resRascon2020ISA.uDiffNorm)
xlim([0 0.5])
boxplot1.YScale = 'log'

fig1sub2 = subplot(1,2,2)
plot(resFerguson.t,resFerguson.qe_norm, resFerguson2.t,resFerguson2.qe_norm, resCruz2021.t,resCruz2021.qe_norm, resRascon2020IET.t,resRascon2020IET.qe_norm, resRascon2020ISA.t,resRascon2020ISA.qe_norm)
fig1sub2.YScale = 'log'
legend('Ferguson 2025 Tuning 1','Ferguson 2025 Tuning 2','Cruz 2021','Rascon 2020 IET','Rascon 2020 ISA','Interpreter','Latex','Orientation','vertical')
xlabel('time (s)')
ylabel('$$||q(t) - q_d(t)||$$','Interpreter','Latex')

% create smaller axes in top right, and plot on it
boxplot2 = axes('Position',[.62 .17 .1 .2])
box on
plot(resFerguson.t,resFerguson.qe_norm, resFerguson2.t,resFerguson2.qe_norm, resCruz2021.t,resCruz2021.qe_norm, resRascon2020IET.t,resRascon2020IET.qe_norm, resRascon2020ISA.t,resRascon2020ISA.qe_norm)
xlim([0 0.5])
boxplot2.YScale = 'log'


%% Functions
function [value, isTerminal, direction] = stopEvent(t,x,stop_func,ph_func)
    % Stop integration when zero detected    
    isTerminal = 1;
    
    q = x(1:2);
    xp = x(5:6);
    phi = x(7);
    ph = ph_func(q,xp,phi);
    value = stop_func(q,ph,phi);
    
    % detect all zeros
    direction = 0;
end

function A = constructAmatrix(q,p,N,dMdqi)
    A = zeros(N);
    for i=1:N
        A(:,i) = feval(dMdqi{i},q)*p;
    end
end

function Sb = constructSbMatrix(q,p,N,S)
    Sb = zeros(N);
    for i=1:N
        ei = zeros(N,1);
        ei(i) = 1;
        Sb(:,i) = S(q,ei)*p;
    end
end

% define saturation functions
function out = sat_vec(x,p,delta)
    out = zeros(size(x));
    for i=1:length(x)
        out(i) = sat(x(i),p,delta);
    end
end

function out = sat(x,p,delta)
    if norm(x) < delta
        out = sign(x)*(abs(x)^p);
    else
        out = delta^p*sign(x);
    end
end