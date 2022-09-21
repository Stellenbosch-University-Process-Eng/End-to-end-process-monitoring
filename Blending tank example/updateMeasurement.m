function [y, sensor_fault] = updateMeasurement(t, x, v, y, sp, d, p, measurement, sensor_fault)
    for i = 1: length(measurement.fields)
        f = measurement.fields{i};
        current = measurement.(f).function(t, x, sp, d, v) + measurement.(f).noise_var*randn;
        
        if strcmp(sensor_fault.(f).state, 'Stuck')
            current = y.(f).Data(end); 

        elseif strcmp(sensor_fault.(f).state, 'Drift')
            current = current + sensor_fault.(f).drift;
            sensor_fault.(f).drift = sensor_fault.(f).drift + sensor_fault.(f).drift_rate;

        elseif strcmp(sensor_fault.(f).state, 'Bias')
            current = current + sensor_fault.(f).bias;

        end

        y.(f).Time = [y.(f).Time t];
        y.(f).Data = [y.(f).Data current];
    end
end

