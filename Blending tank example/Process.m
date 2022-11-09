function x = Process(x, u, d, f, t)

   % Using hardcoded RK4
    x0 = struct2vec(x, t.i);  % Initial values as a vector
    k1 = ODEs(t.time(t.i),          x0,             u, d, f, x.parameters);
    k2 = ODEs(t.time(t.i) + t.dt/2, x0 + t.dt/2*k1, u, d, f, x.parameters);
    k3 = ODEs(t.time(t.i) + t.dt/2, x0 + t.dt/2*k2, u, d, f, x.parameters);
    k4 = ODEs(t.time(t.i) + t.dt,   x0 + t.dt*k3,   u, d, f, x.parameters);
    xvec = (x0 + t.dt/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45

    % Generate the necessary structures
    x = vec2struct(xvec, x, u, d, f, t.time(t.i+1), x.parameters, t.i+1);
end

function dxdt = ODEs(t, xvec, u, d, f, p)
    x = vec2struct(xvec, p.x_empty, u, d, f, t, p, 1);
    
    ddt.m = d.C0(t)*x.F0 - x.C*x.F;
    ddt.V = x.F0 + x.FW - x.F;
    
    ddt.xv = (u.xv - x.xv)/p.tau;
    
    % Check if valve is stuck, which is only an issue if controlled
    % If the valve is forced during startup / shutdown, a bypass can be used
    if (strcmp(f.valveFW.state, 'Stuck')) && (u.control)
        % ddt.xv = 0; % Stuck
        ddt.xv = 1; % Fail open
    
    % Maintain fraction valve opening in [0, 1]
    elseif (x.xv == 0) && (ddt.xv < 0)
        ddt.xv = 0;
    elseif (x.xv == 1) && (ddt.xv > 0)
        ddt.xv = 0;
    end

    ddt.parameters.fields = p.fields;
    dxdt = struct2vec(ddt, 1);
end

function x = vec2struct(xvec, x, u, d, f, t, p, idx)
    for i = 1:length(p.fields)
        x.(p.fields{i})(idx) = xvec(i);
    end
    x.xv(idx) = min(max(x.xv(idx), 0), 1); % ~, valve fraction opening, limit between 0 and 1

    % Add intermediate variables to structure
    x.C(idx)  = x.m(idx)/x.V(idx);
    x.L(idx)  = x.V(idx)/p.A;
    x.FW(idx) = p.cv*x.xv(idx);
    x.F0(idx) = u.F0*d.F0(t);
    x.F(idx)  = u.F*p.kv*real(sqrt(x.L(idx)));
end

function xvec = struct2vec(x, idx)
% Will always take the last element in the x.(field) array 
% i.e., x.(field)(end)
    xvec = zeros(length(x.parameters.fields), 1);
    for i = 1:length(x.parameters.fields)
        xvec(i) = x.(x.parameters.fields{i})(idx);
    end
end



