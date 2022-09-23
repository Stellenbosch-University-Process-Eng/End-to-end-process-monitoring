function [r, t] = SupervisoryControl(r, m, y, t)
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
        r.setpoints.C = [r.setpoints.C nan];
    
        % Special actions associated with this regime
        t.time(end) = t.time(end) + r.Shutdown.period;
            
        % Switching to next regime
        r.regime = 'PrepStartup';
        
    elseif strcmp(r.regime, 'PrepStartup')
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 0;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 0;
        
        % Set-points
        r.setpoints.C = [r.setpoints.C nan];
    
        % Special actions associated with this regime
        % Remove all fault flags from components after shut has been completed
        for i = 1:length(r.components.fields)
            cf = r.components.fields{i};
            if r.components.(cf).faultFlag
                r.components.(cf).faultFlag = false; % Remove maintenance flag after fixing error
                r.components.(cf).commision = t.time(end);  % Signal to fault modules the commission time of a component
            end
        end
        
        % Switching to next regime
        r.Startup.time = [r.Startup.time t.time(end)];
        r.regime = 'Startup';
        
    elseif strcmp(r.regime, 'Startup')
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 1;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 0;
        
        % Set-points
        r.setpoints.C = [r.setpoints.C nan];
    
        % Special actions associated with this regime
        % None
        
        % Switching to next regime
        if y.L.data(end) >= r.Startup.levelThreshold
            r.regime = 'Running';
        end
        
    elseif strcmp(r.regime, 'Running')
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = -1;
        r.components.valveF0.position = 1;
        r.components.valveF.position  = 1;
        
        % Set-points
        r.setpoints.C = [r.setpoints.C 0.3];
    
        % Special actions associated with this regime
        % None for now. Could include flagging components for next maintenance
        
        % Switching to next regime
        if (t.time(end) - r.Startup.time(end)) > r.Running.plannedMaintenancePeriod
            % Planned maintenance
            r.regime = 'Shutdown';
        
        elseif (t.time(end) - r.Startup.time(end) > 3600) && (m.components.C.alarm(end) == 1)
            % Currently, if the sensor is flagged as faulty, immediately shut
            % down plant and perform maintenance
            
            r.components.C.faultFlag = true;
            r.regime = 'Shutdown';
        end
    
    elseif strcmp(r.regime, 'Shutdown')
        % Valve positions (0 = fully closed, 1 = fully open, -1 = controlled)
        r.components.valveFW.position = 0;
        r.components.valveF0.position = 0;
        r.components.valveF.position  = 1;
        
        % Set-points
        r.setpoints.C = [r.setpoints.C nan];
    
        % Special actions associated with this regime
        % None for now. Could include flagging components for next maintenance
        
        % Switching to next regime
        if y.L.data(end) <= r.Shutdown.levelThreshold
            r.regime = 'Shut';
        end
        
    else
        error('No such regime');
    
    end
end