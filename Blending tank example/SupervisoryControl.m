function r = SupervisoryControl(r, m, y, t)
    % r consists of 
    %  Regime: specifies the current operating regime (Shut, Startup, Running, or Shutdown)
    %  components: provides details on sensors, actuators and other components
    %  Setpoints: provides set-points for controllers
    %  Shut / Startup / Running / Shutdown: provides parameters and logs
    %  specifying the operation of the different operating regimes.
    
    if strcmp(r.regime, 'Shut')

        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 0;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 0;
        
        % Set-points
        r.setpoints.C(t.i) = nan;
    
        % Special actions associated with this regime
        % None for now

        % Switching to next regime
        r.regime = 'Startup';
        r.Startup.time(end+1) = t.time(t.i);
        
    elseif strcmp(r.regime, 'Startup')
        
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 1;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 0;
        
        % Set-points
        r.setpoints.C(t.i) = nan;
    
        % Special actions associated with this regime
        % None
        
        % Switching to next regime
        if y.L.data(t.i) >= r.Startup.levelThreshold
            r.regime = 'Running';
        end
        
    elseif strcmp(r.regime, 'Running')
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = -1;
        r.components.valveF0.position = 1;
        r.components.valveF.position  = 1;
        
        % Set-points
        r.setpoints.C(t.i) = r.Running.Csp;
    
        % Special actions associated with this regime
        % None for now.
        
        % Switching to next regime
        if t.time(t.i) > (r.PlannedShuts + 1) * r.PlannedMaintenancePeriod % Planned maintenance
            r.regime = 'Shutdown';
            r.ShutType = r.MaintenanceCycle{mod(r.PlannedShuts, length(r.MaintenanceCycle)) + 1};
            r.PlannedShuts = r.PlannedShuts + 1;

        elseif (t.time(t.i) - r.Startup.time(end) > 3600) % Check for component alarms
            for i = 1:length(r.components.fields)
                cf = r.components.fields{i}; % Current component field
                if (m.components.(cf).alarm(t.i) == 1)  % If any component sounds an alarm...
                    r.components.(cf).faultFlag = true; % ...mark that component as faulty...
                    r.regime = 'Shutdown';              % ... and initiate unplanned maintenance
                    r.ShutType = 'Unplanned';
                end    
            end
            
        end
    
    elseif strcmp(r.regime, 'Shutdown')
        
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 0;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 1;
        
        % Set-points
        r.setpoints.C(t.i) = nan;
    
        % Special actions associated with this regime: 
        % None for now

        % Switching to next regime
        if y.L.data(t.i) <= r.Shutdown.levelThreshold
            r.regime = 'Shut';
        end
        
    else
        error('No such regime');
    
    end
end