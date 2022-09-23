function f = Fault(f, x, r, t)
    % SENSOR FAULTS
    
    % Trigger the developed fault if t > fault time and it hasn't been
    % triggered yet
    if (t.time(end) > f.fault.time) && (~f.fault.triggered)
        f.C.state = 'Drift'; %'None', 'Stuck', 'Drift', 'Bias'
        f.fault.triggered = true;

    elseif strcmp(r.regime, 'Shut') && (r.components.C.faultFlag)
        f.C.state = 'None';    % The supervisory control flagged "C" to be fixed
    end

    % Adjust drift for each sensor which is currently drifting
    for i = 1: length(f.fields)
        if strcmp(f.(f.fields{i}).state, 'Drift')
            f.(f.fields{i}).drift = f.(f.fields{i}).drift + f.(f.fields{i}).driftRate*t.dt;
        end
    end

    % PROCESS FAULTS
    f.valveFW.state = 'None';
end