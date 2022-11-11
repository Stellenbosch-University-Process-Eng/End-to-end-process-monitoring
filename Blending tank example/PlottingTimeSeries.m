%% Plot results
% Prepare functions / variables specific for plotting
clc
LoadCase = 3;
filename_version = '-v5';
load(['Case',num2str(LoadCase), filename_version])

figure(7)


C1 = {'#d7191c','#2c7bb6'};
C2 = {'#1b9e77','#d95f02','#7570b3'};

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

PlotShut.Shade.t = t.time(~isnan(t.time))/3600/24;
PlotShut.Shade.y = (regime(~isnan(t.time)) ~= 3);
PlotShut.idx = find(regime == 1);
PlotShut.Text.t = ( t.time(PlotShut.idx) - r.MinimumShutDownTime )/3600/24;
for i = 1:length(PlotShut.idx)
    PlotShut.Text.Type{i} = shut_type{PlotShut.idx(i)}(1);
end

PlotFault.Shade.t  = t.time(~isnan(t.time))/3600/24;
PlotFault.Shade.C  = faulty_sensor(~isnan(t.time));
PlotFault.Shade.FW = faulty_valve(~isnan(t.time));
PlotFault.Color.C = C1{2};
PlotFault.Color.FW = C2{2};

PlotFault.idx.C = find( (faulty_sensor(2:end) - faulty_sensor(1:end-1)) > 0);
PlotFault.Text.C.t = t.time(PlotFault.idx.C)/3600/24;

PlotFault.idx.FW = find( (faulty_valve(2:end) - faulty_valve(1:end-1)) > 0);
PlotFault.Text.FW.t = t.time(PlotFault.idx.FW)/3600/24;

PlotFault.Alarm.t = t.time(~isnan(t.time))/3600/24;
PlotFault.Alarm.C = m.components.C.alarm(~isnan(t.time));
PlotFault.Alarm.FW = m.components.valveFW.alarm(~isnan(t.time));

% ============
clf
colororder(C1)
cl = tiledlayout('flow','TileSpacing','none');

nexttile(cl)
skip = 25;
line1 = plot(t.time(1:skip:end)/3600/24, x.C(1:skip:end), '^', 'MarkerSize', 3,...
             'DisplayName','$C$');
hold on
line2 = plot(y.C.time(1:skip:end)/3600/24, y.C.data(1:skip:end), 'o', 'MarkerSize', 3,...
             'DisplayName','$C_m$');
line9 = plot(t.time/3600/24, r.setpoints.C,'k--', 'LineWidth', 1.5,...
             'DisplayName','$C^{sp}$');
axis([0 t.tmax/3600/24 -0.2 0.8])
a = gca(); a.YTick = 0:0.2:0.8; a.XTick = [];

[line3, line4, line5, line6] = FaultShading(PlotFault);
ShutDownShading(PlotShut);

nexttile(cl)
colororder(C2)
yyaxis left
line7 = plot(t.time/3600/24, x.L, 'LineWidth', 1.5,...
             'DisplayName','$L$');
xlabel('Time (days)','Interpreter','latex')
axis([0 t.tmax/3600/24 0 6])
a = gca(); a.YTick = 0 : 0.8 : 2.4; a.YColor = 'k';

yyaxis right
line8 = plot(t.time/3600/24, x.xv, 'LineWidth', 1.5,...
             'DisplayName','$x_v$');
axis([0 t.tmax/3600/24 -1.0 1])
a = gca(); a.YTick = 0 : 0.25 : 1; a.YColor = 'k';

%FaultShading(PlotFault);
ShutDownShading(PlotShut, 0);

%hL = legend([line1,line2,line7,line8,line3,line4,line5,line6],'NumColumns',4); 
hL = legend([line9,line1,line2,line7,line3,line4,line8,line5,line6],'NumColumns',3); 
hL.Layout.Tile = 'South'; hL.FontSize = 11;
% legend('True concentration', 'Measured concentration',...
%        '','Sensor fault', 'Sensor alarm',...
%        'Valve fault', 'Valve alarm',...
%        'Location','southoutside','NumColumns', 3)

set(figure(7),'WindowStyle','normal')
set(figure(7), 'Units', 'centimeters');
set(figure(7), 'position', [2 2 8.4 10]);

save(['Case',num2str(r.Case), filename_version]) 
saveas(gcf, ['../Figures/Process-Data-Case',num2str(r.Case), filename_version,'.fig'])
saveas(gcf, ['../Figures/Process-Data-Case',num2str(r.Case), filename_version,'-sized.png'])
disp('Complete')

