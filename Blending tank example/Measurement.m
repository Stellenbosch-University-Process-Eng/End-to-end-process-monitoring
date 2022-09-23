function y = Measurement(y, x, d, fs, t)
    for i = 1: length(y.fields)
        cf = y.fields{i};
        current = y.(cf).function(t.time(end), x, d) + y.(cf).noiseVar*randn;
        
        if strcmp(fs.(cf).state, 'Stuck')
            current = y.(cf).data(end); 

        elseif strcmp(fs.(cf).state, 'Drift')
            current = current + fs.(cf).drift;
            
        elseif strcmp(fs.(cf).state, 'Bias')
            current = current + fs.(cf).bias;

        end

        y.(cf).time = [y.(cf).time t.time(end)];
        y.(cf).data = [y.(cf).data current];
    end
end

