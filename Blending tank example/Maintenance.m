function [r, f, t] = Maintenance(r, f, t)
    % Perform maintenance for unplanned / planned shut downs
    ShutDownTime = r.MinimumShutDownTime;

    % Cycle over all components to identify components requiring
    % replacement and to determine the total shutdown time
    for i = 1:length(r.components.fields)
        cf = r.components.fields{i}; % Current component field
            
        if r.components.(cf).faultFlag ...            % Check any flagged components, or
        || strcmp(r.ShutType, r.components.(cf).type) % check components matching current planned maintenance type
            
            % Add the time taken to check component to shutdown time
            ShutDownTime = ShutDownTime + r.components.(cf).CheckComponentTime;

            % Replace faulty components (regardless of flagged)
            if ~strcmp(f.(cf).state, 'None')
                f.(cf).state = 'None'; % Set the fault state to "None"
                f.(cf).drift = 0; % Reset drift (if any)
                f.(cf).RunTime = 0; % Reset the amount of time component has run
                r.components.(cf).faultFlag = false; % remove the fault flag
                ShutDownTime = ShutDownTime + r.components.(cf).ReplaceComponentTime;
            end
        end
    end
    
    % Move time forward for the duration of the shutdown (proportional to number of inspections)
    t.time(t.i+1) = min(t.time(t.i) + ShutDownTime, t.tmax); % Limit to maximum timepoint
    
end