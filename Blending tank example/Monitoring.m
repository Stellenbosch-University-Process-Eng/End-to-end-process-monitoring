function m = Monitoring(m, y, t)

    if (m.Training) && (t.Time(end) >= m.trainingTime)
        X = [];
        for i = 1:length(m.y_fields)
            f = m.y_fields{i};
            X = [X y.(f).Data(110:end)']; % Remove the first part of transient data
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

        m.Component.C.warning = zeros(1, length(m.T));
        m.Component.C.alarm = false;

        m.Training = false;
        
    else
        X = [];
        for i = 1:length(m.y_fields)
            f = m.y_fields{i};
            X = [X y.(f).Data(end)]; 
        end
        X = (X - m.model.mX)./m.model.sX;
        
        % Monitoring statistics
        T = X*m.model.Q;
        T2 = diag(T * m.model.iSig * T');
        SPE = diag((X - T*m.model.Q') * (X - T*m.model.Q')');
        
        if (T2 > m.hyperparam.T2_threshold) || (SPE > m.hyperparam.SPE_threshold)
            m.Component.C.warning = [m.Component.C.warning 1];
        else
            m.Component.C.warning = [m.Component.C.warning 0];
        end

        if (sum(m.Component.C.warning(end-59:end)) / 60) > 0.8
            m.Component.C.alarm = true;
        end

        m.statistic.T = [m.statistic.T T];
        m.statistic.T2 = [m.statistic.T2 T2];
        m.statistic.SPE = [m.statistic.SPE SPE];

    end


end
