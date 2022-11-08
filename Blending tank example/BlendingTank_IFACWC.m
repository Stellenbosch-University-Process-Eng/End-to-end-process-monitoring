%% Initialize
clc
clear
clf
rng(2)
tic
%% Time span (t)
t.dt = 1;       % s
t.tmax = 4*3600;  % s, simulation time, 8*3600 works well for now

% Pre-allocate for speed
N = t.tmax / t.dt + 1;  % Max array size, if no shutdowns occur
t.time = NaN(N, 1); % Create column array of NaN values

t.i = 1;    % Current time index
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
% Component specific parameters
r.components.fields = {'valveF0','valveFW','valveF', 'C',     'C0',    'F0',    'FW',    'F',     'L'};
component_types     = {'Valve',  'Valve',  'Valve',  'Sensor','Sensor','Sensor','Sensor','Sensor','Sensor'};
for i = 1:length(r.components.fields)
    cf = r.components.fields{i}; % Current component field
    r.components.(cf).type = component_types{i};
    r.components.(cf).faultFlag = false; % Whether the component has been flagged as faulty
    r.components.(cf).CheckComponentTime = 15*60;   % s, time taken to check component
    r.components.(cf).ReplaceComponentTime = 30*60; % s, time taken to replace component
end
clear component_types

% Maintenance parameters
r.MinimumShutDownTime = 600;   % s, minimum shutdown period
r.PlannedMaintenancePeriod = 3*3600;  % s, time before planned maintenance
r.MaintenanceCycle = {'Valve', 'Sensor'}; % Cycle for planned maintenance actions
r.PlannedShuts = 0; % Number of planned shuts that have occured (used to estimate position in cycle)

% Regime specific parameters
r.Shutdown.levelThreshold = 0.001;   % m, level at which to switch from "Shutdown" to "Shut"
r.Startup.levelThreshold = 0.5;     % m, level at which to switch from "Startup" to "Running"
r.Startup.time = 0;
r.Running.Csp = 0.3;    % Concentration set-point during the "Running" regime

% Initialize supervisory control
r.setpoints.C = NaN(N,1);
r.regime = 'Startup'; % Operating regime when simulation starts

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
    
    % Pre-allocate for speed
    x.(x.parameters.fields{i}) = NaN(N,1);
end
for i = 1:length(x.parameters.intFields)
    x.parameters.x_empty.(x.parameters.intFields{i}) = [];

    % Pre-allocate for speed
    x.(x.parameters.intFields{i}) = NaN(N,1);
end


%% Faults (f)
% CDF for bathtub failure model: alpha ~ minimum probability for failure, L = max lifetime
BathtubCDF = @(t, alpha, L) (16*(1-alpha)*(t/L-1/2).^5 + alpha*(t/L-1/2) + 1/2);
f.fields = r.components.fields;

% Each component has an associated fault state, 
% which may be "None", "Bias", "Drift" or "Stuck" for sensors,
% and "None" or "Stuck" for valves
% Each sensor also has a drift rate associated with the "Drift" fault, 
% such that the drift at time "t" is given by drift = drift_rate * (t - t(incipient fault))
% The bias simply gives the amount by which the sensor is offset for the "Bias" fault

% Faults for concentration measurement
f.C.F = @(t) BathtubCDF(t, 0.05, 1.6e4); % CDF of failure rate; alpha ~ minimum probability for failure, L = max lifetime        
f.C.fault_type = 'Drift';   % If a fault occurs, it will be a drift fault
f.C.drift = 0;
f.C.driftRate = 0.0001;

% All other sensor faults; none will be introduced for this example
f.C0.F = @(t) 0; f.C0.fault_type = 'None';
f.F0.F = @(t) 0; f.F0.fault_type = 'None';
f.FW.F = @(t) 0; f.FW.fault_type = 'None';
f.F.F = @(t) 0;  f.F.fault_type = 'None';
f.L.F = @(t) 0;  f.L.fault_type = 'None';

% Valve faults
f.valveFW.F = @(t) BathtubCDF(t, 0.1, 2*3600); 
f.valveFW.fault_type = 'Stuck';

% All other valve faults; none will be introduced for this example
f.valveF0.F = @(t) 0; f.valveF0.fault_type = 'None';
f.valveF.F = @(t) 0;  f.valveF.fault_type = 'None';

% Determine hazard function for each component, set initial fault state and run time
for i = 1: length(f.fields)
    cf = f.fields{i};   % Current component field
    % Hazard function for failure, see https://en.wikipedia.org/wiki/Failure_rate#Failure_rate_in_the_discrete_sense
    % The component fails if a uniform random number between 0 and 1 is smaller than f.C.hazard(t)
    f.(cf).hazard = @(time) ( f.(cf).F(time + t.dt) - f.(cf).F(time) ) ./ ( ( 1-f.(cf).F(time) ) * t.dt );

    % Set all component fault states and run times equal to zero
    f.(cf).state = 'None';
    f.(cf).RunTime = 0; 

end

%% Measurement (y)
% List of measurement variables
% Each measurement has an associated function which is used to calculate
% the measurement from the process variables, as well as a noise variance
% which specifies the variance of the normally distruted noise added to
% each measurement
y.fields = {'C','C0','F0','FW', 'F', 'L'};

% Concentration in the tank
y.C.function = @(t, x, d, idx) x.C(idx);
y.C.noiseVar = 0.01;

% Inlet concentration
y.C0.function = @(t, x, d, idx) d.C0(t);
y.C0.noiseVar = 0.01;

% Inlet flowrate
y.F0.function = @(t, x, d, idx) x.F0(idx);
y.F0.noiseVar = 0.002;

% Water flowrate
y.FW.function = @(t, x, d, idx) x.FW(idx);
y.FW.noiseVar = 0.002;

% Liquid level
y.L.function = @(t, x, d, idx) x.L(idx);
y.L.noiseVar = 0.002;

% Outlet flowrate
y.F.function = @(t, x, d, idx) x.F(idx);
y.F.noiseVar = 0.002;

% Initialize measurements (pre-allocate for speed)
for i = 1:length(y.fields)
    y.(y.fields{i}).time = NaN(N,1);
    y.(y.fields{i}).data = NaN(N,1);
end


%% Monitoring (m)
% m.yFields are the measurements which are used during monitoring, to
% create a warning or an alarm
% m.components.fields are the different components that may have an alarm
% associated with it. For example, even though the component "valveFW"
% isn't a measured variable, a specific combination of measurements may
% indicate that the feed water valve is malfunctioning, which will yield an
% alarm in m.components.valveFW

m.yFields = {'C','C0','F0','FW', 'F', 'L'};
m.components.fields = r.components.fields;

% Model hyperparameters
m.hyperparam.nComponents = 2;
m.hyperparam.T2_threshold = 30;
m.hyperparam.SPE_threshold = 20;

% Current alarms or warnings on any component
% Alarms are passed to the supervisory control layer
% Pre-allocate for speed
for i = 1:length(m.components.fields)
    cf = m.components.fields{i};
    m.components.(cf).alarm = NaN(N,1);
    m.components.(cf).warning = NaN(N,1);
end

% Specify model training time
m.training = true;      % Determines if monitoring method is still training
m.trainingTime = 2000;  % Time taken to train monitoring method

% Pre-allocate for speed
m.statistic.T = NaN(N, m.hyperparam.nComponents);
m.statistic.T2 = NaN(N,1);
m.statistic.SPE = NaN(N,1);

%% Economic model
% Pre-allocate for speed
econ.KPI.values = NaN(N,1);
econ.KPI.function = @(r, x, idx) exp( -40*(x.C(idx) - r.setpoints.C(idx-1)).^2 );

%% Initial conditions
% Initialize the process state variables
t.time(t.i) = 0;
x.m(t.i) = 0.5; % kg, initial solute concentration in tank
x.V(t.i) = 0.25;   % m3, initial liquid volume in tank
x.xv(t.i) = 0.5; % ~, initial fraction valve opening
x.v(t.i) = 0;   % 1/s, initial valve velocity

%% Simulate
while t.time(t.i) < t.tmax
    t.time(t.i + 1) = t.time(t.i) + t.dt;
    
    r = SupervisoryControl(r, m, y, t);    
    if strcmp(r.regime, 'Shut')
        [r, f, t] = Maintenance(r, f, t);  % Maintenance can update "t" during "Shut"
    end
    u = RegulatoryControl(u, y, r, t);
    x = Process(x, u, d, f, t);
    f = Fault(f, x, t);
    y = Measurement(y, x, d, f, t);
    m = Monitoring(m, y, r, t);
    econ = Economic(econ, r, x, t);
    disp(t.time(t.i)/t.tmax)
    
    t.i = t.i + 1;
end
disp('Done')



%% Plot results
subplot(2,2,1)
plot(y.C.time/3600, y.C.data, '.', ...
     t.time/3600, x.C, '.', ...
     t.time/3600, r.setpoints.C,'k--', ...
     t.time(m.components.C.alarm == 1)/3600, 0.35*ones(sum(m.components.C.alarm), 1),'r|',...
     'LineWidth', 2)
xlabel('Time (s)'); ylabel('Concentration'); 
legend('Measured C', 'Actual C', 'Set-point C','Location','best')
axis([0 t.tmax/3600 0 0.8])

subplot(2,2,2)
plot(t.time/3600, x.L, '.')
xlabel('Time'); ylabel('Liquid level');
% plot(t.Time, econ.KPI.Values, 'LineWidth', 2)
% xlabel('Time'); ylabel('KPI');

subplot(2,2,3)
plot(y.C0.time/3600, y.C0.data,'.', ...
     t.time/3600, d.C0(t.time),'.')
xlabel('Time'); ylabel('C0');

subplot(2,2,4)
plot(t.time/3600,    x.F0, ...
     t.time/3600,    x.FW, ...
     t.time/3600,    x.F, 'LineWidth', 2)
xlabel('Time'); ylabel('F');
legend('F0', 'FW', 'F');
toc