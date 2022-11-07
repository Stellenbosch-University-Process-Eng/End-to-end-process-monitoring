function x = Process(x, u, d, f, t)

   % Using hardcoded RK4
    x0 = struct2vec(x);          % Initial values as a vector
    k1 = ODEs(t.time(end-1),          x0,             u, d, f, x.parameters);
    k2 = ODEs(t.time(end-1) + t.dt/2, x0 + t.dt/2*k1, u, d, f, x.parameters);
    k3 = ODEs(t.time(end-1) + t.dt/2, x0 + t.dt/2*k2, u, d, f, x.parameters);
    k4 = ODEs(t.time(end) ,           x0 + t.dt*k3,   u, d, f, x.parameters);
    xvec = (x0 + t.dt/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45

    % Generate the necessary structures
    x = vec2struct(xvec(end,:), x, u, d, f, t.time(end), x.parameters);
end

function dxdt = ODEs(t, xvec, u, d, f, p)
    x = vec2struct(xvec, p.x_empty, u, d, f, t, p);
    
    ddt.m = d.C0(t)*x.F0 - x.C*x.F;
    ddt.V = x.F0 + x.FW - x.F;
    ddt.xv = x.v;
    ddt.v = (1/p.tau^2)*(u.xv - x.xv) - 2*p.xi/p.tau*x.v;
    
    % Maintain fraction valve opening in [0, 1]
    if (x.xv == 0) && (ddt.v < 0)
        ddt.xv = 0;
        ddt.v = 0;
    elseif (x.xv == 1) && (ddt.v > 0)
        ddt.xv = 0;
        ddt.v = 0;
    end
    
    ddt.parameters.fields = p.fields;
    dxdt = struct2vec(ddt);
end

function x = vec2struct(xvec, x, u, d, f, t, p)
    for i = 1:length(p.fields)
        x.(p.fields{i}) = [x.(p.fields{i}); xvec(i)];
    end
    x.xv(end) = min(max(x.xv(end), 0), 1); % ~, valve fraction opening, limit between 0 and 1

    % Add intermediate variables to structure
    x.C  = [x.C  x.m(end)/x.V(end)];
    x.L  = [x.L  x.V(end)/p.A];
    x.FW = [x.FW p.cv*x.xv(end)];
    x.F0 = [x.F0 u.F0*d.F0(t)];
    x.F  = [x.F  u.F*p.kv*real(sqrt(x.L(end)))];
end

function xvec = struct2vec(x)
% Will always take the last element in the x.(field) array 
% i.e., x.(field)(end)
    xvec = zeros(length(x.parameters.fields), 1);
    for i = 1:length(x.parameters.fields)
        xvec(i) = x.(x.parameters.fields{i})(end);
    end
end



