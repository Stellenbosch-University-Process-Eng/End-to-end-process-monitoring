function y = Measurement(y, x, d, f, t)
    for i = 1: length(y.fields)
        cf = y.fields{i};
        current = y.(cf).function(t.time(end), x, d) + y.(cf).noiseVar*randn;
        
        if strcmp(f.(cf).state, 'Stuck')
            current = y.(cf).data(end); 

        elseif strcmp(f.(cf).state, 'Drift')
            current = current + f.(cf).drift;
            
        elseif strcmp(f.(cf).state, 'Bias')
            current = current + f.(cf).bias;

        end

        y.(cf).time(end+1) = t.time(end);
        y.(cf).data(end+1) = current;
    end
end

