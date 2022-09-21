function y = initMeasurement(x0, sp, d, p, measurement)
    % Populate measurements with zeros. The zero measurements will be used
    % for the first "intermediateVariables" call
    for i = 1:length(measurement.fields)
        f = measurement.fields{i};
        y.(f).Time = 0;
        y.(f).Data = 0;
    end
    
    v = intermediateVariables(0, x0, y, sp, d, p);
    for i = 1:length(measurement.fields)
        f = measurement.fields{i};
        y.(f).Data = measurement.(f).function(0, x0, sp, d, v);
    end
end

