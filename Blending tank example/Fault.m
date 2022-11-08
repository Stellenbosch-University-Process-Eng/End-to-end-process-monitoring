function f = Fault(f, x, t)
    
    for i = 1: length(f.fields)
        cf = f.fields{i};   % Current component field
        
        % Increase the component's running time
        f.(cf).RunTime = f.(cf).RunTime + t.dt;

        % Trigger faults according to the probability for failure defined
        % by hazard function
        if rand < f.(cf).hazard(f.(cf).RunTime)
           f.(cf).state = f.(cf).fault_type; 
        end
        
        % Adjust drift for each sensor which is currently drifting
        if strcmp(f.(cf).state, 'Drift')
            f.(cf).drift = f.(cf).drift + f.(cf).driftRate*t.dt;
        end
    end

end