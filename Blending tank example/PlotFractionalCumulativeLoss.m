%% Compare cumulative KPI
C2 = {'#1b9e77','#d95f02','#7570b3'};
colororder(C2);
figure(6)
clf
filename_version = '-v5';
for j = 1:3
    load(['Case',num2str(j), filename_version]) 
    cKPI{j} = cumsum(econ.KPI.values,'omitnan');
end
colororder({'#1b9e77','#d95f02','#7570b3'});

subplot(2,1,1)
plot(t.time/3600/24, cKPI{1},'--',...
     t.time/3600/24, cKPI{2},...
     t.time/3600/24, cKPI{3},'LineWidth', 1.5);

ylabel('Cumulative profit')
xlabel('Time (days)')
legend('No intervention', ...
       'Planned maintenance only', ...
       'Planned and unplanned maintenance', ...
       'Location','northwest')

subplot(2,1,2)
plot(t.time/3600/24, 0*t.time + 1,'--',...
     t.time/3600/24, cKPI{2}./cKPI{1},...
     t.time/3600/24, cKPI{3}./cKPI{1}, 'LineWidth', 1.5);

ylabel('Fractional cumulative loss/gain')
xlabel('Time (days)')
legend('No intervention', ...
       'Planned maintenance only', ...
       'Planned and unplanned maintenance', ...
       'Location','northwest')

axis([0 t.tmax/3600/24 0.8 1.6])
%a = gca(); a.YTick = a.YTick(a.YTick <= 1);

saveas(gcf, ['../Figures/Cumulative-Fractional-Loss-Gain', filename_version,'.fig'])
saveas(gcf, ['../Figures/Cumulative-Fractional-Loss-Gain', filename_version,'.tif'])
