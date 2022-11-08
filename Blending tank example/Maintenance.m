function [r, f, t] = Maintenance(r, f, t)
    % Perform maintenance for unplanned / planned shut downs
    ShutDownTime = 0;

    % Cycle over all components
    for i = 1:length(r.components.fields)
        cf = r.components.fields{i}; % Current component field
            
        if r.components.(cf).faultFlag ...            % Check any flagged components, or
        || strcmp(r.ShutType, r.components.(cf).type) % check components matching current planned maintenance type
            
            % Add the time taken to check component to shutdown time
            ShutDownTime = ShutDownTime + r.components.(cf).

            % if those components are also faulty (regardless of flagged)
            ~strcmp(f.(cf).state, 'None'
            f.(cf).state = 'None'; % Set the fault state to zero (replace component)
            f.(cf).drift = 0; % Reset drift (if any)
            r.components.(cf).faultFlag = false; % remove the fault flag
            r.components.(cf).commision = t.time(end); % s, time when new component was commissioned
        end

    end
    
    % Move time forward for the duration of the shutdown 
    % (proportional to number of inspections)
    t.time(end) = t.time(end) + ShutDownPeriod;
    

end