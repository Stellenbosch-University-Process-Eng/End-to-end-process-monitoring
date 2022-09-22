function u = Control(u, y, r, t)
u.F0 = r.Components.valveF.position;
u.F = r.Components.valveF.position;

if r.Components.valveFW.position == -1   % Apply control
    error = r.Setpoint.C - y.C.Data(end);
    u.intError = u.intError + error*t.dt;
    u.xv = -u.PI.K*(error + u.intError/u.PI.tauI);
    
else
    u.xv = r.Components.valveFW.position;    % Will be either 0 (closed) or 1 (open)
    u.intError = 0;

end