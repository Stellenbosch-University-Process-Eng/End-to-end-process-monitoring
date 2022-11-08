function y = Measurement(y, x, d, f, t)
    for i = 1: length(y.fields)
        cf = y.fields{i};
        current = y.(cf).function(t.time(t.i+1), x, d, t.i+1) + y.(cf).noiseVar*randn;
        
        if strcmp(f.(cf).state, 'Stuck')
            current = y.(cf).data(end); 

        elseif strcmp(f.(cf).state, 'Drift')
            current = current + f.(cf).drift;
            
        elseif strcmp(f.(cf).state, 'Bias')
            current = current + f.(cf).bias;

        end

        y.(cf).time(t.i+1) = t.time(t.i+1);
        y.(cf).data(t.i+1) = current;
    end
end

