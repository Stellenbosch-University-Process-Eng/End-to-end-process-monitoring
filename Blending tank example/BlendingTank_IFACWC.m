%% Initialize
clc
clear
clf
rng(1)

%% Time span (t)
t.dt = 1;       % s
t.tmax = 6*3600;  % s, simulation time

%% Disturbance variables (d)
% Create stochastic inlet flowrate and concentrations over time
F0 = 0; C0 = 0;
tspan = 0: t.dt : t.tmax;
for i = 2:length(tspan)
    F0(i) = 0.99*F0(i-1) + 0.00015*randn;
    C0(i) = 0.999*C0(i-1) + 0.005*randn;
end
F0 = F0 + 0.01; 
C0 = C0 + 1;    

d.F0 = griddedInterpolant(tspan, F0);
d.C0 = griddedInterpolant(tspan, C0);
clear F0 C0 tspan

%% Supervisory control (r)
r.components.fields = {'valveF0','valveFW','valveF', 'C','C0','F0','FW', 'F', 'L'};
for i = 1:length(r.components.fields)
    r.components.(r.components.fields{i}).faultFlag = false;
    r.components.(r.components.fields{i}).commision = 0;
end

r.Shutdown.period = 3600;   % s, length of a shut down
r.Shutdown.levelThreshold = 0.001;   % m, level at which to switch from "Shutdown" to "Shut"
r.Startup.levelThreshold = 0.5;     % m, level at which to switch from "Startup" to "Running"
r.Running.plannedMaintenancePeriod = 2*3600;  % s, time before planned maintenance
r.Startup.time = 0;

r.regime = 'PrepStartup';
r.setpoints.C = nan;

%% Regulatory control (u)
u.PI.K = 0.1;      % m3/kg, controller gain
u.PI.tauI = 10;    % s, controller time constant

%% Process (x)
% Define process parameters
x.parameters.A  = 0.5;     % m2, mixing tank cross-sectional area
x.parameters.tau = 10;     % s, valve time constant
x.parameters.xi = 5;       % ~, valve damping coefficient
x.parameters.cv = 0.1;     % m3/s, control valve coefficient
x.parameters.kv = 0.06;    % m2.5/s, drainage valve coefficient

% List of state variables. Create an empty array for each state value,
% which will be used in the Simulate function
x.parameters.fields = {'m', 'V', 'xv', 'v'};            % Fields for state variables
x.parameters.intFields = {'C', 'L','FW', 'F0', 'F'};    % Fields for intermediate variables

% Create a structure with all variable fields empty, useful when calling ODEs
for i = 1:length(x.parameters.fields)
    x.parameters.x_empty.(x.parameters.fields{i}) = [];
end
for i = 1:length(x.parameters.intFields)
    x.parameters.x_empty.(x.parameters.intFields{i}) = [];
end


%% Faults (f)
% SENSOR FAULTS
% Each measurement has an associated fault state, 
% which may be "None", "Bias", "Drift" or "Stuck". 
% Each sensor also has a drift rate associated with the "Drift" fault, 
% such that the drift at time "t" is given by drift = drift_rate * (t - t(incipient fault))
% The bias simply gives the amount by which the sensor is offset for the
% "Bias" fault

f.fields = {'C','C0','F0','FW', 'F', 'L'};
% Faults for concentration measurement
f.C.state = 'None';
f.C.drift = 0;
f.C.driftRate = 0.0001;
f.C.bias = 0.2;

% All other faults; none will be introduced for this example
f.C0.state = 'None';
f.F0.state = 'None';
f.FW.state = 'None';
f.F.state = 'None';
f.L.state = 'None';
f.KPI.state = 'None';

% Parameters specifying when a fault might occur
f.fault.time = 2185;
f.fault.triggered = false;

% PROCESS FAULTS
% Initialize process fault
f.valveFW.state = 'None';

%% Measurement (y)
% List of measurement variables
% Each measurement has an associated function which is used to calculate
% the measurement from the prcoess variables, as well as a noise variance
% which specifies the variance of the normally distruted noise added to
% each measurement
y.fields = {'C','C0','F0','FW', 'F', 'L'};

% Concentration in the tank
y.C.function = @(t, x, d) x.C(end);
y.C.noiseVar = 0.01;

% Inlet concentration
y.C0.function = @(t, x, d) d.C0(t);
y.C0.noiseVar = 0.01;

% Inlet flowrate
y.F0.function = @(t, x, d) x.F0(end);
y.F0.noiseVar = 0.002;

% Water flowrate
y.FW.function = @(t, x, d) x.FW(end);
y.FW.noiseVar = 0.002;

% Liquid level
y.L.function = @(t, x, d) x.L(end);
y.L.noiseVar = 0.002;

% Outlet flowrate
y.F.function = @(t, x, d) x.F(end);
y.F.noiseVar = 0.002;

% Initialize measurements
for i = 1:length(y.fields)
    y.(y.fields{i}).time = [];
    y.(y.fields{i}).data = [];
end


%% Monitoring (m)
% Measurements to include in monitoring model
m.yFields = {'C','C0','F0','FW', 'F', 'L'};

% Model hyperparameters
m.hyperparam.nComponents = 2;
m.hyperparam.T2_threshold = 30;
m.hyperparam.SPE_threshold = 20;

% Specify model training time
m.training = true;      % Determines if monitoring method is still trainign
m.trainingTime = 2000;  % Time taken to train monitoring method

% Current flags on any component. 
% Alarms are passed to the supervisory control layer
m.component.C.alarm = false;

%% Economic model
econ.KPI.values = nan;
econ.KPI.function = @(r, x) exp( -40*(x.C(end) - r.setpoints.C(end)).^2 );

%% Simulate
% Initialize the process state variables
x.m = 0.5; % kg, initial solute concentration in tank
x.V = 0.25;   % m3, initial liquid volume in tank
x.xv = 0.5; % ~, initial fraction valve opening
x.v = 0;   % 1/s, initial valve velocity
for i = 1:length(x.parameters.intFields)
    x.(x.parameters.intFields{i}) = nan;
end

t.time = 0;
while t.time(end) < t.tmax
    t.time = [t.time t.time(end)+t.dt];
    
    [r, t] = SupervisoryControl(r, m, y, t);    % Supervisory control can update "t" during "Shut"
    u = RegulatoryControl(u, y, r, t);
    x = Process(x, u, d, f, t);
    f = Fault(f, x, r, t);
    y = Measurement(y, x, d, f, t);
    m = Monitoring(m, y, t);
    econ = Economic(econ, r, x);

    disp(t.time(end)/t.tmax)
end
disp('Done')

%% Plot results
subplot(2,2,1)
plot(y.C.time, y.C.data, '.', ...
     t.time, x.C, '.', ...
     t.time, r.setpoints.C,'k--', ...
     t.time(m.components.C.alarm == 1), 0.35*ones(sum(m.components.C.alarm), 1),'r|',...
     'LineWidth', 2)
xlabel('Time (s)'); ylabel('Concentration'); 
legend('Measured C', 'Actual C', 'Set-point C','Location','best')
axis([0 t.tmax 0 0.8])

subplot(2,2,2)
plot(t.time, x.L, '.')
xlabel('Time'); ylabel('Liquid level');
% plot(t.Time, econ.KPI.Values, 'LineWidth', 2)
% xlabel('Time'); ylabel('KPI');

subplot(2,2,3)
plot(y.C0.time, y.C0.data,'.', ...
     t.time, d.C0(t.time),'.')
xlabel('Time'); ylabel('C0');

subplot(2,2,4)
plot(t.time,    x.F0, ...
     t.time,    x.FW, ...
     t.time,    x.F, 'LineWidth', 2)
xlabel('Time'); ylabel('F');
legend('F0', 'FW', 'F');
