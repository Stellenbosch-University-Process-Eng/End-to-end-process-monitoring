function fs = SensorFault(fs, x, r, t)
    
    % Trigger the developed fault if t > fault time and it hasn't been
    % triggered yet
    if (t.time(end) > fs.fault.time) && (~fs.fault.triggered)
        fs.C.state = 'Drift'; %'None', 'Stuck', 'Drift', 'Bias'
        fs.fault.triggered = true;

    elseif strcmp(r.regime, 'Shut') && (r.components.C.faultFlag)
        fs.C.state = 'None';    % The supervisory control flagged "C" to be fixed
    end

    % Adjust drift for each sensor which is currently drifting
    for i = 1: length(fs.fields)
        f = fs.fields{i};
        if strcmp(fs.(f).state, 'Drift')
            fs.(f).drift = fs.(f).drift + fs.(f).driftRate*t.dt;
        end
    end
end