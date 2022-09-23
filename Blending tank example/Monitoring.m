function m = Monitoring(m, y, t)

    if (m.training) && (t.time(end) >= m.trainingTime)
        X = [];
        for i = 1:length(m.yFields)
            X = [X y.(m.yFields{i}).data(110:end)']; % Remove the first part of transient data
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

        m.components.C.warning = zeros(1, length(m.statistic.T));
        m.components.C.alarm = zeros(1, length(m.statistic.T));

        m.training = false;
        
    elseif ~m.training
        X = [];
        for i = 1:length(m.yFields)
            X = [X y.(m.yFields{i}).data(end)]; 
        end
        X = (X - m.model.mX)./m.model.sX;
        
        % Monitoring statistics
        T = X*m.model.Q;
        T2 = diag(T * m.model.iSig * T');
        SPE = diag((X - T*m.model.Q') * (X - T*m.model.Q')');
        
        if (T2 > m.hyperparam.T2_threshold) || (SPE > m.hyperparam.SPE_threshold)
            m.components.C.warning = [m.components.C.warning 1];
        else
            m.components.C.warning = [m.components.C.warning 0];
        end

        if (sum(m.components.C.warning(end-59:end)) / 60) > 0.8
            m.components.C.alarm = [m.components.C.alarm 1];
        else
            m.components.C.alarm = [m.components.C.alarm 0];
        end

        m.statistic.T = [m.statistic.T; T];
        m.statistic.T2 = [m.statistic.T2; T2];
        m.statistic.SPE = [m.statistic.SPE; SPE];

    end


end
