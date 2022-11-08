function econ = Economic(econ, r, x, t)


econ.KPI.values(t.i + 1) = econ.KPI.function(r, x, t.i+1);
