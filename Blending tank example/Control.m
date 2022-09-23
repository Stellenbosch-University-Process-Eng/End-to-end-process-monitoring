function u = Control(u, y, r, t)
u.F0 = r.components.valveF0.position;
u.F = r.components.valveF.position;

if r.components.valveFW.position == -1   % Apply control
    error = r.setpoints.C(end) - y.C.data(end);
    u.intError = u.intError + error*t.dt;
    u.xv = -u.PI.K*(error + u.intError/u.PI.tauI);
    
else   % Will be either 0 (closed) or 1 (open)
    u.xv = r.components.valveFW.position;
    u.intError = 0;

end