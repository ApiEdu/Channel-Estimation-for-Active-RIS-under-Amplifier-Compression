function mdpi_style(figname, widthCM)
%MDPI_STYLE  Apply MDPI figure requirements and export the current figure.
%
%   mdpi_style('fig_nmse_floor_ibo')        % single-column figure, 12 cm
%   mdpi_style('fig_scaling', 17)           % full-width two-panel figure
%
%   MDPI asks for at least 300 dpi for colour figures. Telecom is a
%   single-column layout, so text that looked fine in the IEEE two-column
%   version comes out small here; the font sizes below compensate.
%
%   Call this as the LAST line of each figN script, after the plot is drawn.

if nargin < 2 || isempty(widthCM), widthCM = 12; end   % cm
heightCM = widthCM * 0.68;

fig = gcf;
set(fig, 'Units','centimeters', 'Position',[2 2 widthCM heightCM], ...
         'PaperUnits','centimeters', 'PaperSize',[widthCM heightCM], ...
         'PaperPosition',[0 0 widthCM heightCM], 'Color','w');

for ax = findall(fig,'Type','axes').'
    set(ax, 'FontSize',11, 'FontName','Helvetica', 'LineWidth',0.9, ...
            'Box','on', 'TickDir','in', 'XGrid','on','YGrid','on', ...
            'GridLineStyle',':', 'GridAlpha',0.35);
    set(get(ax,'XLabel'), 'FontSize',12);
    set(get(ax,'YLabel'), 'FontSize',12);
    set(get(ax,'Title'),  'FontSize',12, 'FontWeight','normal');
    lg = get(ax,'Legend'); if ~isempty(lg), set(lg,'FontSize',10); end
end
set(findall(fig,'Type','line'), 'LineWidth',1.8);

% 300 dpi PNG. exportgraphics needs R2020a+; print() is the fallback.
try
    exportgraphics(fig, [figname '.png'], 'Resolution',300, ...
                   'BackgroundColor','white');
catch
    print(fig, [figname '.png'], '-dpng', '-r300');
end
fprintf('exported %s.png at 300 dpi (%g x %.1f cm)\n', figname, widthCM, heightCM);
end
