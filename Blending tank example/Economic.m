function econ = Economic(econ, r, x)


econ.KPI.values = [econ.KPI.values econ.KPI.function(r, x)];
