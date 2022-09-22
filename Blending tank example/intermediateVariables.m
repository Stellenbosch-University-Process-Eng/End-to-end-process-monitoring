function x = intermediateVariables(x, u, d, fp, p)
    x.C = x.m/x.V;
    x.L = x.V/p.A;
    x.FW = p.cv*x.xv;
    x.F  = p.kv*sqrt(x.L);
end
