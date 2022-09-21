function [T, T2, SPE, faulty, monitor] = initMonitoring(y, monitor)
    X = [y.C.Data' y.C0.Data' y.F0.Data' y.FW.Data' y.F.Data' y.L.Data'];
    X = X(110:end,:);   % Remove the first part of transient data
    monitor.mX = mean(X);
    monitor.sX = std(X);
    X = (X - monitor.mX)./monitor.sX;
    
    [monitor.Q, Sig] = eigs((X' * X)/length(X), monitor.nComponents);
    monitor.iSig = inv(Sig);
    
    T = X*monitor.Q;

    T2 = diag(T * monitor.iSig * T');
    SPE = diag((X - T*monitor.Q') * (X - T*monitor.Q')');
    
    faulty = zeros(1, length(T));
end