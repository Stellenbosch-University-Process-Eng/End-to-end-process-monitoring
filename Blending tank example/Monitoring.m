function [T, T2, SPE, faulty] = Monitoring(y, monitor)
    X = [y.C.Data(end) y.C0.Data(end) y.F0.Data(end) y.FW.Data(end) y.F.Data(end) y.L.Data(end)];
    X = (X - monitor.mX)./monitor.sX;
    T = X*monitor.Q;
    
    T2 = T * monitor.iSig * T';
    SPE = (X - T*monitor.Q') * (X - T*monitor.Q')';
    
    if (T2 > monitor.T2_threshold) || (SPE > monitor.SPE_threshold)
        faulty = 1;
    else
        faulty = 0;
    end
end