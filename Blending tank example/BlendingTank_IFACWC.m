%% Initialize
clc
clear
clf
rng(1)

%% Time span
t = 0:7200;    % s, over one day

%% Parameters
% Define process parameters
p.A  = 0.5;     % m2, mixing tank cross-sectional area
p.tauI = 1;     % s, valve time constant
p.xi = 5;       % ~, valve damping coefficient
p.K = 0.1;      % m3/kg, controller gain
p.tauI = 10;    % s, controller time constant
p.cv = 0.1;     % m3/s, control valve coefficient
p.kv = 0.06;    % m2.5/s, drainage valve coefficient

% List of state variables. Create an empty array for each state value,
% which will be used in the Simulate function
p.state_fields = {'m', 'V', 'xv', 'v'};
for i = 1:length(p.state_fields)
    p.x_empty.(p.state_fields{i}) = [];
end
%% Measurement parameters
% List of measurement variables
% Each measurement has an associated function which is used to calculate
% the measurement from the prcoess variables, as well as a noise variance
% which specifies the variance of the normally distruted noise added to
% each measurement
y.fields = {'C','C0','F0','FW', 'F', 'L', 'KPI'};

% Concentration in the tank
y.C.function = @(t, x, sp, d, v) x.C(end);
y.C.noise_var = 0.01;

% Inlet concentration
y.C0.function = @(t, x, sp, d, v) d.C0(t);
y.C0.noise_var = 0.01;

% Inlet flowrate
y.F0.function = @(t, x, sp, d, v) d.F0(t);
y.F0.noise_var = 0.002;

% Water flowrate
y.FW.function = @(t, x, sp, d, v) v.FW(end);
y.FW.noise_var = 0.002;

% Liquid level
y.L.function = @(t, x, sp, d, v) x.L(end);
y.L.noise_var = 0.002;

% Outlet flowrate
y.F.function = @(t, x, sp, d, v) v.F(end);
y.F.noise_var = 0.002;

% Key performance indicator
y.KPI.function = @(t, x, sp, d, v) exp( -40*(x.C(end) - sp.C(t)).^2 );
y.KPI.noise_var = 0;
    
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

%% Monitoring parameters
m.nComponents = 2;
m.T2_threshold = 30;
m.SPE_threshold = 20;

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

%% Set-point changes
% Specify set-point changes
r.C = @(t) 0.5 + 0.*(t>(3600*0.5));

%% Initialize
% Initialize the process variables
x.m = 0.5; % kg, initial solute concentration in tank
x.V = 1;   % m3, initial liquid volume in tank
x.xv = 0.5; % ~, initial fraction valve opening
x.v = 0;   % 1/s, initial valve velocity
x = intermediateVariables(x, u, d, fp, p);

%% Integrate ODEs
% Initialize the simulation
%y = initMeasurement(x, sp, d, p, y); % Initialize measurements
shut = [];  % s, vector containing time values where plant was shut down
fault_time = 2185;

disp('Simulation started')

i = 1;
while t(i) < t(end)
    i = i+1;
    x = Simulate(x, u, d, fp, t(i-1:i), p);
    [y, fs] = updateMeasurement(t(i), x, v, y, sp, d, p, y, fs);
    
    % Monitoring related activities
    if t(i) == 2000 % Train monitoring model
        [T, T2, SPE, faulty, m] = initMonitoring(y, m);
        
    elseif t(i) > 2000 % Continue with monitoring
        [T(i,:), T2(i), SPE(i), faulty(i)] = Monitoring(y, m);

        if (sum(faulty(i-59:i)) / 60) > 0.8
            disp('Shut down');
            shut = [shut t(i)];
            fs.C.state = 'None';
        end
    end
    
    % Fault related activities
    % Trigger a sensor fault
    if (t(i) > fault_time) && isempty(shut)
        fs.C.state = 'Drift'; %'None', 'Stuck', 'Drift', 'Bias'
    end
    disp(i/length(t))
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
