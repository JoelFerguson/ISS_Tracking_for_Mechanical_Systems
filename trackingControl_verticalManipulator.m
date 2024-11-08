%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ISS tracking control for mechanical systems
% Description: Simulation files for the paper 'Input-to-State Stable 
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

%% Simulation settings
sim.isSave = 0;        % 0 = Dont save results; 1 = Save results
sim.disturbance = 1;   % 0 = No disturbance; 1 = Disturbance active

% Simulation initial conditions
sim.q0 = [0 0].';   % configuration (q)
sim.p0 = [-1 2].';  % momentum (p)
sim.ph0 = [0 0].';  % momentum estimate (\hat p)
sim.t_end = 20;      % simulation length

% Observer parameters
obs.kappa = 10;

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
[t,x] = ode45(sim.ode,[0 sim.t_end],sim.x0,sim.options);
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
res.t = t;
res.q = x(:,1:sys.N);
res.p0 = x(:,sys.N+1:2*sys.N);
res.xp = x(:,2*sys.N+1:3*sys.N);
res.phi = x(:,3*sys.N+1);

% Preallocate space for results
res.H = zeros(length(t),1);
res.Ho = zeros(length(t),1);
res.He = zeros(length(t),1);
res.W = zeros(length(t),1);
res.p0h = zeros(length(t),2);
res.p = zeros(length(t),2);
res.ph = zeros(length(t),2);
res.eig_Sb = zeros(length(t),1);
res.switchCond = zeros(length(t),1);
res.qe = zeros(length(t),2);
res.pe = zeros(length(t),2);
res.qd = zeros(length(t),2);
res.qe_norm = zeros(length(t),1);

% Evaluate functions of interest for each time step
for i=1:length(t)
    % Unpack state for the ith time step
    t = res.t(i);
    q = res.q(i,:).';
    p0 = res.p0(i,:).';
    xp = res.xp(i,:).';
    phi = res.phi(i);

    % Compute coordinate transforms and error coordinates
    p = obs.Ti(q)\p0;
    ph = obs.ph(q,xp,phi);

    % Evaluate functions of interest
    res.H(i) = sys.H(q,p0);
    res.p(i,:) = p.';
    res.ph(i,:) = ph.';
    res.p0h(i,:) = (obs.Ti(q)*ph).';
    res.Ho(i) = obs.Hp(p,ph);
    res.u(i,:) = sim.u(q,ph,phi,t).';
    res.qe(i,:) = ctrl.qe(q,t).';
    res.pe(i,:) = ctrl.pe(q,ph,t).';
    res.He(i) = ctrl.He(res.qe(i,:).', res.pe(i,:).');
    res.W(i) = res.Ho(i) + res.He(i);
    res.qd(i,:) = ctrl.qd(t).';
    res.qe_norm(i) = norm(res.qe(i,:));

    % compute max eigenvalues of symm[Sb]
    res.eig_Sb(i) = max(eig(0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).')));
    res.switchCond(i) = obs.switchCond(q,ph,phi);
end

% tracking error plots
fig1 = figure("Name","Momentum estimate");
plot(res.t,res.p,res.t,res.ph,'--')
legend('$$p_1$$','$$p_2$$','$$\hat{p}_1$$','$$\hat{p}_2$$','Interpreter','Latex','NumColumns',4)
xlabel('time (s)')
ylabel('Momentum estimate')

fig2 = figure("Name","Performance");
subplot(2,1,1)
plot(res.t,log(res.W))
xlabel('time (s)')
ylabel('Joint storage function')
legend("$log(W)$")

subplot(2,1,2)
plot(res.t,res.phi)
xlabel('time (s)')
ylabel('$\phi$')
legend("Observer gain")

fig3 = figure("Name","Tracking error");
plot(res.t,res.q,res.t,res.qd,'--')
legend('$$q_1$$','$$q_2$$','$$q_{d1}$$','$$q_{d2}$$','Interpreter','Latex','NumColumns',4)
xlabel('time (s)')
ylabel('Configuration')

%% Save figure
if sim.isSave
    results_name = sprintf("observerResults_kappa_%i_disturbance_%i",obs.kappa*100,sim.disturbance);
    save(results_name,'res');
end

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