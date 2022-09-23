function y = Measurement(y, x, d, fs, t)
    for i = 1: length(y.fields)
        f = y.fields{i};
        current = y.(f).function(t.time(end), x, d) + y.(f).noiseVar*randn;
        
        if strcmp(fs.(f).state, 'Stuck')
            current = y.(f).data(end); 

        elseif strcmp(fs.(f).state, 'Drift')
            current = current + fs.(f).drift;
            
        elseif strcmp(fs.(f).state, 'Bias')
            current = current + fs.(f).bias;

        end

        y.(f).time = [y.(f).time t.time(end)];
        y.(f).data = [y.(f).data current];
    end
end

