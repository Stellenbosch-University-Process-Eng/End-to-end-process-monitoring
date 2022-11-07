%%
clc
clear
clf

Bathtub = @(t, alpha, L) (16*(1-alpha)*(t/L-1/2).^5 + alpha*(t/L-1/2) + 1/2);
alpha = 0.05; L = 1000;

F = @(t) Bathtub(t, alpha, L);
t = linspace(0, 1000, 1000);
dt = t(2);

h = @(t) ( F(t + dt) - F(t) ) ./ ( (1-F(t))*dt );
%plot(t, h(t))

N = 10000;
Survive = zeros(N, length(t));
Survive(:, 1) = 1;
Failure = t(end) * ones(N,1);

for i = 2:length(t)
    rnd = rand(N,1);
    Survive(:,i) = Survive(:,i-1) .* ( rnd > h(t(i)) );
    Failure( Survive(:,i) ~= Survive(:,i-1)) = t(i);
end

histogram(Failure,'Normalization','cdf')
hold on
plot(t, F(t), "LineWidth", 2);
hold off