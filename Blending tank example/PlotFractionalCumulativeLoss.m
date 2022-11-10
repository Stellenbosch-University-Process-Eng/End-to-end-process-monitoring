%% Compare cumulative KPI
figure(4)
clf
filename_version = '-v4';
for j = 1:3
    load(['Case',num2str(j), filename_version]) 
    cKPI{j} = cumsum(econ.KPI.values,'omitnan');
end
colororder({'#1b9e77','#d95f02','#7570b3'});

plot(t.time/3600/24, cKPI{2}./cKPI{1},...
     t.time/3600/24, cKPI{3}./cKPI{1}, 'LineWidth', 1.5);

ylabel('Fractional cumulative loss/gain')
xlabel('Time (days)')
legend('Planned maintenance only', ...
       'Planned and unplanned maintenance', ...
       'Location','northwest')

axis([0 t.tmax/3600/24 0.95 1.2])
%a = gca(); a.YTick = a.YTick(a.YTick <= 1);

saveas(gcf, ['Cumulative-Fractional-Loss-Gain', filename_version,'.fig'])
saveas(gcf, ['Cumulative-Fractional-Loss-Gain', filename_version,'.tif'])
