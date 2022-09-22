function fs = SensorFault(fs, x, r, t)
    
    % Trigger the developed fault if t > fault time and it hasn't been
    % triggered yet
    if (t.Time(end) > fs.fault.time) && (~fs.fault.triggered)
        fs.C.state = 'Drift'; %'None', 'Stuck', 'Drift', 'Bias'
    elseif (fs.fault.triggered) && (r.replace.C)
        fs.C.State = 'None';    % The supervisory control flagged "C" to be fixed
    end

    % Adjust drift for each sensor which is currently drifting
    for i = 1: length(fs.fields)
        f = fs.fields{i};
        if strcmp(fs.(f).state, 'Drift')
            fs.(f).drift = fs.(f).drift + fs.(f).drift_rate*t.dt;
        end
    end
end