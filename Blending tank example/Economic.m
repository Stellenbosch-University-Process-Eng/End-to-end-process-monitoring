function econ = Economic(econ, r, x, t)



if strcmp(r.regime, 'Shut')
    if strcmp(r.ShutType, 'Unplanned')
        econ.KPI.values(t.i + 1) = -250*r.ShutDownTime/3600;
    else
        econ.KPI.values(t.i + 1) = -100*r.ShutDownTime/3600;
    end
else
    econ.KPI.values(t.i + 1) = econ.KPI.function(r, x, t.i+1);
end
