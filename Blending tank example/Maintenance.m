function [r, f] = Maintenance(r, f, t)
    % Perform maintenance for unplanned / planned shut downs
    
    % Cycle over all components
    for i = 1:length(r.components.fields)
        cf = r.components.fields{i}; % Current component field
            
        if r.components.(cf).faultFlag ... % Replace any flagged components, or
        || (     strcmp(r.ShutType, r.components.(cf).type) ... % Replace components matching current planned maintenance type ...
             && ~strcmp(f.(cf).state, 'None') )                 % if those components are also faulty (regardless of flagged)

            f.(cf).state = 'None'; % Set the fault state to zero (replace component)
            f.(cf).drift; % Reset drift (if any)
            r.components.(cf).faultFlag = false; % remove the fault flag
            r.components.(cf).commision = t.time(end); % s, time when new component was commissioned
        end

    end

end