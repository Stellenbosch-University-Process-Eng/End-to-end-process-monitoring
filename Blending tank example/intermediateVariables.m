function v = intermediateVariables(t, x, y, sp, d, p)
    v.FW = p.cv*x.x;
    v.F  = p.kv*sqrt(x.L);

    v.e = sp.C(t) - y.C.Data(end);
    
    v.xsp = -p.K*(v.e + x.I/p.tauI);
end
