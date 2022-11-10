function FaultShading(PlotFault)

    % Indicate shutdowns
    hold on
    a = gca(); h = a.YLim(2) - a.YLim(1); b = a.YLim(1);
    shade = fill([PlotFault.Shade.t' PlotFault.Shade.t(end) 0], 0.10*h*[PlotFault.Shade.C' 0 0]+ b,'k');
    shade.FaceAlpha = 0.2;
    shade.FaceColor = PlotFault.Color.C;
    shade.EdgeColor = PlotFault.Color.C;
    shade.LineStyle = '-';

    text(PlotFault.Text.C.t, 0.10*h + b + zeros(length(PlotFault.Text.C.t),1), ' C',...
        'VerticalAlignment','top')
    
    plot(PlotFault.Alarm.t, 0.3*h*PlotFault.Alarm.C + b - 0.2*h, '|', ...
        'Color', PlotFault.Color.C, 'LineWidth', 5);

    shade = fill([PlotFault.Shade.t' PlotFault.Shade.t(end) 0], 0.20*h*[PlotFault.Shade.FW' 0 0]+ b,'k');
    shade.FaceAlpha = 0.2;
    shade.FaceColor = PlotFault.Color.FW;
    shade.EdgeColor = PlotFault.Color.FW;
    shade.LineStyle = '-';
    
    text(PlotFault.Text.FW.t, 0.20*h + b + zeros(length(PlotFault.Text.FW.t),1), ' FW',...
        'VerticalAlignment','top')
    
    plot(PlotFault.Alarm.t, 0.4*h*PlotFault.Alarm.FW + b - 0.2*h, '|', ...
        'Color', PlotFault.Color.FW, 'LineWidth', 5);


    hold off

end