function y = Measurement(y, x, fs, t)
    for i = 1: length(y.fields)
        f = y.fields{i};
        current = y.(f).function(t.Time(end), x, sp, d) + y.(f).noise_var*randn;
        
        if strcmp(fs.(f).state, 'Stuck')
            current = y.(f).Data(end); 

        elseif strcmp(fs.(f).state, 'Drift')
            current = current + fs.(f).drift;
            
        elseif strcmp(fs.(f).state, 'Bias')
            current = current + fs.(f).bias;

        end

        y.(f).Time = [y.(f).Time t.Time(end)];
        y.(f).Data = [y.(f).Data current];
    end
end

