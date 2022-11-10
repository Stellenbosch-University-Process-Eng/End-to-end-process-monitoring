function ShutDownShading(PlotShut)

    % Indicate shutdowns
    hold on
    a = gca(); h = a.YLim(2) - a.YLim(1); b = a.YLim(1);
    shade = fill([PlotShut.Shade.t' PlotShut.Shade.t(end) 0], h*[PlotShut.Shade.y' 0 0]+b, 'k');
    shade.FaceAlpha = 0.1;
    shade.EdgeColor = 'None';
    text(PlotShut.Text.t, 0.95*h + b + zeros(length(PlotShut.Text.t),1), PlotShut.Text.Type)
    hold off

end