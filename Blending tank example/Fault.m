function f = Fault(f, x, t)
    
    for i = 1: length(f.fields)
        cf = f.fields{i};   % Current component field
        
        % Trigger faults according to the probability for failure defined
        % by hazard function
        if rand < f.(cf).hazard(t.time(end))
           f.(cf).state = f.(cf).fault_type; 
        end
        
        % Adjust drift for each sensor which is currently drifting
        if strcmp(f.(cf).state, 'Drift')
            f.(cf).drift = f.(cf).drift + f.(cf).driftRate*t.dt;
        end
    end

    % PROCESS FAULTS
    f.valveFW.state = 'None';
end