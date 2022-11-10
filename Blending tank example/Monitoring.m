function m = Monitoring(m, y, r, t)

    if (m.training) && (t.time(t.i+1) >= m.trainingTime)
        % Runs once, when training time is over

        X = [];
        for i = 1:length(m.yFields)
            % X = [X y.(m.yFields{i}).data(110:end)']; % Remove the first part of transient data
            X(:,end+1) = y.(m.yFields{i}).data(110:t.i+1)'; % Remove the first part of transient data
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
            m.components.(cf).warning(1:t.i+1) = 0;
            m.components.(cf).alarm(1:t.i+1) = 0;
        end
        m.training = false;
        
    elseif ~(m.training) && (strcmp(r.regime, 'Running'))
        % Set all alarms / warning to zero at first
        for i = 1:length(m.components.fields)
            cf = m.components.fields{i};
            m.components.(cf).warning(t.i+1) = 0;
            m.components.(cf).alarm(t.i+1) = 0;
        end

        if m.MonitoringActive
        
            % Add the latest measurement
            X = NaN(1, length(m.yFields));
            for i = 1:length(m.yFields)
                X(i) = y.(m.yFields{i}).data(t.i+1); 
            end
            X = (X - m.model.mX)./m.model.sX;
            
            % Monitoring statistics
            T = X*m.model.Q;
            T2 = diag(T * m.model.iSig * T');
            SPE = diag((X - T*m.model.Q') * (X - T*m.model.Q')');
            
            m.statistic.T(t.i+1,:) = T;
            m.statistic.T2(t.i+1) = T2;
            m.statistic.SPE(t.i+1) = SPE;
    
            % Check if the composition alarm tripped
            if (T2 > m.hyperparam.T2_threshold) || (SPE > m.hyperparam.SPE_threshold)
                if abs(y.C.data(t.i+1) - m.model.mX(1)) < 4*m.model.sX(1) 
                    % If there is a problem, but the sensor reading is close to
                    % SP, then assume sensor is faulty...
                    m.components.C.warning(t.i+1) = 1;
                else
                    % ... else assume valve is faulty
                    m.components.valveFW.warning(t.i+1) = 1;
                end
            end
    
            % ### Note the hard-coded values
            alarm_w = 36;  % Monitoring window size
            alarm_f = 0.8; % Fraction of warnings in window to sound alarm
            if (sum(m.components.C.warning(t.i - alarm_w+2: t.i+1)) / alarm_w) > alarm_f
                m.components.C.alarm(t.i+1) = 1;
            end
            if (sum(m.components.valveFW.warning(t.i - alarm_w+2: t.i+1)) / alarm_w) > alarm_f
                m.components.valveFW.alarm(t.i+1) = 1;
            end
        
        else
            if y.C.data(t.i+1) > (m.model.mX(1) - 2*m.model.sX(1))
                m.MonitoringActive = true;
            end
        end
    elseif strcmp(r.regime,'Startup')
        m.MonitoringActive = false;   % Provide time to get to Csp
    
    end


end
