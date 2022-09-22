function x = Process(x, u, d, fp, t)

%     % Using ode45
%     [~, xvec] = ode45(@(t, x) ODEs(t, x, y, sp, d, p), [t.Time(end-1) t.Time(end)], struct2vec(x, p));
    
    % Using hardcoded RK4
    x0 = struct2vec(x);          % Initial values as a vector
    k1 = ODEs(t.Time(end-1),          x0,           u, d, fp, x);
    k2 = ODEs(t.Time(end-1) + t.dt/2, x0 + dt/2*k1, u, d, fp, x);
    k3 = ODEs(t.Time(end-1) + t.dt/2, x0 + dt/2*k2, u, d, fp, x);
    k4 = ODEs(t.Time(end) ,           x0 + dt*k3,   u, d, fp, x);
    xvec = (x0 + t.dt/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45

    % Generate the necessary structures
    x = vec2struct(xvec(end,:), x);
end

function dxdt = ODEs(t, xvec, u, d, fp, x)
    x = vec2struct(xvec, x);
    p = x.parameters;
    
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
    
    ddt.fields = x.fields;
    dxdt = struct2vec(ddt);
end

function x = intermediateVariables(x, u, fp)
    p = x.parameters;
    
    x.C = x.m/x.V;
    x.L = x.V/p.A;
    x.FW = p.cv*x.xv;
    x.F0 = u.F0*d.F0(t);
    x.F  = u.F*p.kv*sqrt(x.L);
end

function x = vec2struct(xvec, x)
    for i = 1:length(x.fields)
        x.(x.fields{i}) = [x.(x.fields{i}); xvec(i)];
    end
    x.xv(end) = min(max(x.xv(end), 0), 1); % ~, valve fraction opening, limit between 0 and 1

    % Add intermediate variables to structure
    x = intermediateVariables(x, u, fp);
end

function xvec = struct2vec(x)
% Will always take the last element in the x.(field) array 
% i.e., x.(field)(end)
    xvec = zeros(length(x.fields), 1);
    for i = 1:length(x.fields)
        xvec(i) = x.(x.fields{i})(end);
    end
end



