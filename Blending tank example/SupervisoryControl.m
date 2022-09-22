function [r, t] = SupervisoryControl(r, m, y, t)
% r consists of 
%  Regime: specifies the current operating regime (Shut, Startup, Running, or Shutdown)
%  Components: provides details on sensors, actuators and other components
%  Setpoints: provides set-points for controllers
%  Shut / Startup / Running / Shutdown: provides parameters and logs
%  specifying the operation of the different operating regimes.

if strcomp(r.Regime, 'Shut')
    % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
    r.Components.valveFW.position = 0;
    r.Components.valveF0.position = 0;
    r.Components.valveF.position  = 0;
    
    % Set-points
    r.Setpoints.C = 0;

    % Special actions associated with this regime
    % Perform maintenance on all flagged components
    for i = 1:length(r.Components.fields)
        f = r.Components.fields{i};
        if r.Components.(f).faultFlag
            r.Components.(f).faultFlag = false; % Remove maintenance flag after fixing error
            r.Components.(f).commision = t.Times(end);  % Signal to fault modules the commission time of a component
        end
    end
    
    % Switching to next regime
    t.Times = [t.Times t.Times(end) + r.Shutdown.Period];
    r.Startup.times = [r.Startup.times t.Times(end)];
    r.Regime = 'Startup';
    
elseif strcomp(r.Regime, 'Startup')
    % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
    r.Components.valveFW.position = 1;
    r.Components.valveF0.position = 0;
    r.Components.valveF.position  = 0;
    
    % Set-points
    r.Setpoints.C = 0;

    % Special actions associated with this regime
    % None
    
    % Switching to next regime
    if y.L >= r.Startup.levelThreshold
        r.Regime = 'Running';
    end
    
elseif strcomp(r.Regime, 'Running')
    % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
    r.Components.valveFW.position = -1;
    r.Components.valveF0.position = 1;
    r.Components.valveF.position  = 1;
    
    % Set-points
    r.Setpoints.C = 0.5;

    % Special actions associated with this regime
    % None for now. Could include flagging components for next maintenance
    
    % Switching to next regime
    if (t - r.Startup.times(end)) > r.Running.plannedMaintenancePeriod
        % Planned maintenance
        r.Regime = 'Shutting';
    
    elseif m.Component.C.faultFlag
        % Currently, if the sensor is flagged as faulty, immediately shut
        % down plant and perform maintenance
        
        r.Component.C.faultFlag = true;
        r.Regime = 'Shutdown';
    end

elseif strcomp(r.regime, 'Shutdown')
    % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
    r.Components.valveFW.position = 0;
    r.Components.valveF0.position = 0;
    r.Components.valveF.position  = 1;
    
    % Set-points
    r.Setpoints.C = 0;

    % Special actions associated with this regime
    % None for now. Could include flagging components for next maintenance
    
    % Switching to next regime
    if y.L <= r.Shutdown.levelThreshold
        r.Regime = 'Shut';
    end
    
else
    error('No such regime');

end