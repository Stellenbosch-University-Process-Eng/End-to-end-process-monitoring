function [x, v] = Simulate(tspan, x, y, sp, d, p)

%     % Using ode45
%     [~, xvec] = ode45(@(t, x) ODEs(t, x, y, sp, d, p), tspan, struct2vec(x, p));
    
    % Using hardcoded RK4
    dt = tspan(end) - tspan(1);     % Time step
    x0 = struct2vec(x, p);          % Initial values as a vector
    k1 = ODEs(tspan(1), x0, y, sp, d, p);
    k2 = ODEs(tspan(1) + dt/2, x0 + dt/2*k1, y, sp, d, p);
    k3 = ODEs(tspan(1) + dt/2, x0 + dt/2*k2, y, sp, d, p);
    k4 = ODEs(tspan(1) + dt, x0 + dt*k3, y, sp, d, p);
    xvec = (x0 + dt/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45

    % Generate the necessary structures
    x = vec2struct(xvec(end,:), x, p);
    v = intermediateVariables(tspan(end), x, y, sp, d, p);

end

function dxdt = ODEs(t, xvec, y, sp, d, p)
    x = vec2struct(xvec, p.x_empty, p);
    v = intermediateVariables(t, x, y, sp, d, p);

    ddt.C = ( d.C0(t)*d.F0(t) - x.C*v.F ) / (p.A*x.L);
    ddt.L = (d.F0(t) + v.FW - v.F) / p.A;
    ddt.x = x.v;
    ddt.v = (p.k*(v.xsp - x.x) - p.s*x.v);
    ddt.I = v.e;
    
    % Maintain fraction valve opening in [0, 1]
    if (x.x == 0) && (ddt.v < 0)
        ddt.x = 0;
        ddt.v = 0;
    elseif (x.x == 1) && (ddt.v > 0)
        ddt.x = 0;
        ddt.v = 0;
    end

    dxdt = struct2vec(ddt, p);
end

function x = vec2struct(xvec, x, p)
    for i = 1:length(p.state_fields)
        x.(p.state_fields{i}) = [x.(p.state_fields{i}); xvec(i)];
    end
    x.x(end) = min(max(x.x(end), 0), 1); % ~, valve fraction opening, limit between 0 and 1
end

function xvec = struct2vec(x, p)
% Will always take the last element in the x.(field) array 
% i.e., x.(field)(end)
    xvec = zeros(length(p.state_fields), 1);
    for i = 1:length(p.state_fields)
        xvec(i) = x.(p.state_fields{i})(end);
    end
end



