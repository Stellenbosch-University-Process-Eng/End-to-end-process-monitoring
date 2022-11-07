function m = Monitoring(m, y, t)

    if (m.training) && (t.time(end) >= m.trainingTime)
        % Runs once, when training time is over

        X = [];
        for i = 1:length(m.yFields)
            % X = [X y.(m.yFields{i}).data(110:end)']; % Remove the first part of transient data
            X(:,end+1) = y.(m.yFields{i}).data(110:end)'; % Remove the first part of transient data
        end
        
        % Center and scale data
        m.model.mX = mean(X);
        m.model.sX = std(X);
        X = (X - m.model.mX)./m.model.sX;
        
        % Train PCA model
        [m.model.Q, Sig] = eigs((X' * X)/length(X), m.hyperparam.nComponents);
        m.model.iSig = inv(Sig);
        
        % Monitoring statistics
        T = X*m.model.Q;
        T2 = diag(T * m.model.iSig * T');
        SPE = diag((X - T*m.model.Q') * (X - T*m.model.Q')');
        
        m.statistic.T = T;
        m.statistic.T2 = T2;
        m.statistic.SPE = SPE;


        % Set current warnings / alarms
        for i = 1:length(m.components.fields)
            cf = m.components.fields{i};
            m.components.(cf).warning = zeros(1, length(m.statistic.T));
            m.components.(cf).alarm = zeros(1, length(m.statistic.T));
        end
        m.training = false;
        
    elseif ~m.training
        
        % Add the latest measurement
        X = [];
        for i = 1:length(m.yFields)
            X(end+1) = y.(m.yFields{i}).data(end); 
        end
        X = (X - m.model.mX)./m.model.sX;
        
        % Monitoring statistics
        T = X*m.model.Q;
        T2 = diag(T * m.model.iSig * T');
        SPE = diag((X - T*m.model.Q') * (X - T*m.model.Q')');
        
        m.statistic.T(end+1) = T;
        m.statistic.T2(end+1) = T2;
        m.statistic.SPE(end+1) = SPE;

        % Set all alarms / warning to zero at first
        for i = 1:length(m.components.fields)
            cf = m.components.fields{i};
            m.components.(cf).warning(end+1) = 0;
            m.components.(cf).alarm(end+1) = 0;
        end
                
        
        % Check if the composition alarm tripped
        if (T2 > m.hyperparam.T2_threshold) || (SPE > m.hyperparam.SPE_threshold)
            m.components.C.warning(end) = 1;
        end

        if (sum(m.components.C.warning(end-59:end)) / 60) > 0.8
            m.components.C.alarm(end) = 1;
        end
        
    end


end
