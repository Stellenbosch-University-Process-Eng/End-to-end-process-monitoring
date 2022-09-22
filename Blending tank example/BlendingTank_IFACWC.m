%% Initialize
clc
clear
clf
rng(1)

%% Time span
t.dt = 1;       % s
t.tmax = 7200;  % s, simulation time

%% Process parameters
% Define process parameters
x.parameters.A  = 0.5;     % m2, mixing tank cross-sectional area
x.parameters.tauI = 1;     % s, valve time constant
x.parameters.xi = 5;       % ~, valve damping coefficient
x.parameters.cv = 0.1;     % m3/s, control valve coefficient
x.parameters.kv = 0.06;    % m2.5/s, drainage valve coefficient

% List of state variables. Create an empty array for each state value,
% which will be used in the Simulate function
x.fields = {'m', 'V', 'xv', 'v'};
for i = 1:length(x.fields)
    x.x_empty.(x.fields{i}) = [];
end
%% Measurement parameters
% List of measurement variables
% Each measurement has an associated function which is used to calculate
% the measurement from the prcoess variables, as well as a noise variance
% which specifies the variance of the normally distruted noise added to
% each measurement
y.fields = {'C','C0','F0','FW', 'F', 'L', 'KPI'};

% Concentration in the tank
y.C.function = @(t, x, sp, d) x.C(end);
y.C.noise_var = 0.01;

% Inlet concentration
y.C0.function = @(t, x, sp, d) d.C0(t);
y.C0.noise_var = 0.01;

% Inlet flowrate
y.F0.function = @(t, x, sp, d) x.F0(end);
y.F0.noise_var = 0.002;

% Water flowrate
y.FW.function = @(t, x, sp, d) x.FW(end);
y.FW.noise_var = 0.002;

% Liquid level
y.L.function = @(t, x, sp, d) x.L(end);
y.L.noise_var = 0.002;

% Outlet flowrate
y.F.function = @(t, x, sp, d) x.F(end);
y.F.noise_var = 0.002;

% Key performance indicator
y.KPI.function = @(t, x, sp, d) exp( -40*(x.C(end) - sp.C(t)).^2 );
y.KPI.noise_var = 0;
    
% Initialize measurements
for i = 1:length(y.fields)
    f = y.fields{i};
    y.(f).Time = [];
    y.(f).Data = [];
end

%% Process faults (fp)
% Initialize process fault
fp.valve.state = 'None';

%% Sensor faults (fs)
% Each measurement has an associated fault state, 
% which may be "None", "Bias", "Drift" or "Stuck". 
% Each sensor also has a drift rate associated with the "Drift" fault, 
% such that the drift at time "t" is given by drift = drift_rate * (t - t(incipient fault))
% The bias simply gives the amount by which the sensor is offset for the
% "Bias" fault

fs.fields = y.fields;

% Faults for concentration measurement
fs.C.state = 'None';
fs.C.drift = 0;
fs.C.drift_rate = 0.0001;
fs.C.bias = 0.2;

% All other faults; none will be introduced for this example
fs.C0.state = 'None';
fs.F0.state = 'None';
fs.FW.state = 'None';
fs.F.state = 'None';
fs.L.state = 'None';
fs.KPI.state = 'None';

% Parameters specifying when a fault might occur
fs.fault.time = 2185;
fs.fault.triggered = false;

%% Monitoring parameters
% Measurements to include in monitoring model
m.y_fields = {'C','C0','F0','FW', 'F', 'L'};

% Model hyperparameters
m.hyperparam.nComponents = 2;
m.hyperparam.T2_threshold = 30;
m.hyperparam.SPE_threshold = 20;

% Specify model training time
m.Training = true;      % Determines if monitoring method is still trainign
m.trainingTime = 2000;  % Time taken to train monitoring method

% Current flags on any component
m.Component.C.faultFlag = false;

%% Disturbance variables

% Create stochastic inlet flowrate and concentrations over time
F0 = 0*t; C0 = 0*t;
for i = 2:length(t)
    F0(i) = 0.99*F0(i-1) + 0.00015*randn;
    C0(i) = 0.999*C0(i-1) + 0.005*randn;
end
F0 = F0 + 0.01; 
C0 = C0 + 1;    

d.F0 = griddedInterpolant(t, F0);
d.C0 = griddedInterpolant(t, C0);
clear F0 C0

%% Supervisory control
r.Components.fields = {'valveF0','valveFW','valveF', 'C','C0','F0','FW', 'F', 'L'};
for i = 1:length(r.Components.fields)
    f = r.Components.fields{i};
    r.Components.(f).faultFlag = false;
    r.Components.(f).commision = 0;
end

r.Shutdown.Period = 3600;   % s, length of a shut down
r.Shutdown.levelThreshold = 0.01;   % m, level at which to switch from "Shutdown" to "Shut"
r.Startup.levelThreshold = 0.5;     % m, level at which to switch from "Startup" to "Running"
r.Running.plannedMaintenancePeriod = 7200;  % s, time before planned maintenance
r.Regime = 'Startup';
r.Startup.times = 0;

%% Regulatory control
u.PI.K = 0.1;      % m3/kg, controller gain
u.PI.tauI = 10;    % s, controller time constant

%% Initialize
% Initialize the process state variables
x.m = 0.5; % kg, initial solute concentration in tank
x.V = 1;   % m3, initial liquid volume in tank
x.xv = 0.5; % ~, initial fraction valve opening
x.v = 0;   % 1/s, initial valve velocity

%% Integrate ODEs
% Initialize the simulation
%y = initMeasurement(x, sp, d, p, y); % Initialize measurements
shut = [];  % s, vector containing time values where plant was shut down

disp('Simulation started')

t.Time = 0;
while t.Time(end) < t.tmax
    t.Time = [t.Time t.Time(end)+t.dt];
    
    [r, t] = SupervisoryControl(r, m, y, t);
    u = Control(u, y, r, t);
    x = Process(x, u, d, fp, t, p);
    fp = ProcessFault(fp, x, t);
    fs = SensorFault(fs, x, t);
    y = Measurement(y, x, fs, t);
    m = Monitoring(m, y, t);
    
    disp(t.Time(end)/t.tmax)
end
disp('Done')

%% Plot results
subplot(2,2,1)
plot(y.C.Time, y.C.Data, '.', ...
     t, x.C, ...
     t, sp.C(t),'k--', ...
     t(faulty == 1), 0.35*ones(sum(faulty), 1),'|',...
     'LineWidth', 2)
xlabel('Time (s)'); ylabel('Concentration'); 
legend('Measured C', 'Actual C', 'Set-point C','Location','best')

subplot(2,2,2)
plot(y.KPI.Time, y.KPI.Data, 'LineWidth', 2)
xlabel('Time'); ylabel('KPI');

subplot(2,2,3)
plot(y.C0.Time, y.C0.Data,'.', ...
     t, d.C0(t), 'LineWidth', 2)
xlabel('Time'); ylabel('C0');

subplot(2,2,4)
plot(y.F0.Time, y.F0.Data,'.', ...
     t, d.F0(t), ...
     y.FW.Time, y.FW.Data,'.', ...
     t, v.FW, 'LineWidth', 2)
xlabel('Time'); ylabel('F');
legend('Measured F0','True F0', 'Measured FW', 'True FW');
