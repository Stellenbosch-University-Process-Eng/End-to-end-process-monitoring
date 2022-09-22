function x = Process(x, u, d, fp, tspan, p)

%     % Using ode45
%     [~, xvec] = ode45(@(t, x) ODEs(t, x, y, sp, d, p), tspan, struct2vec(x, p));
    
    % Using hardcoded RK4
    dt = tspan(end) - tspan(1);     % Time step
    x0 = struct2vec(x, p);          % Initial values as a vector
    k1 = ODEs(tspan(1), x0, u, d, fp, p);
    k2 = ODEs(tspan(1) + dt/2, x0 + dt/2*k1, u, d, fp, p);
    k3 = ODEs(tspan(1) + dt/2, x0 + dt/2*k2, u, d, fp, p);
    k4 = ODEs(tspan(1) + dt, x0 + dt*k3, u, d, fp, p);
    xvec = (x0 + dt/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45

    % Generate the necessary structures
    x = vec2struct(xvec(end,:), x, p);
end

function dxdt = ODEs(t, xvec, u, d, fp, p)
    x = vec2struct(xvec, p.x_empty, p);

    ddt.m = d.C0(t)*d.F0(t) - x.C*x.F;
    ddt.V = d.F0(t) + x.FW - x.F;
    ddt.xv = x.v;
    ddt.v = (1/p.tau^2)*(u.xv - x.xv) - 2*p.xi/p.tau*x.v);
    
    % Maintain fraction valve opening in [0, 1]
    if (x.xv == 0) && (ddt.v < 0)
        ddt.xv = 0;
        ddt.v = 0;
    elseif (x.xv == 1) && (ddt.v > 0)
        ddt.xv = 0;
        ddt.v = 0;
    end

    dxdt = struct2vec(ddt, p);
end

function x = vec2struct(xvec, x, p)
    for i = 1:length(p.state_fields)
        x.(p.state_fields{i}) = [x.(p.state_fields{i}); xvec(i)];
    end
    x.xv(end) = min(max(x.xv(end), 0), 1); % ~, valve fraction opening, limit between 0 and 1

    % Add intermediate variables to structure
    x = intermediateVariables(x, u, d, fp, p);
end

function xvec = struct2vec(x, p)
% Will always take the last element in the x.(field) array 
% i.e., x.(field)(end)
    xvec = zeros(length(p.state_fields), 1);
    for i = 1:length(p.state_fields)
        xvec(i) = x.(p.state_fields{i})(end);
    end
end



