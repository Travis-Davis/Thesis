function figureformat(fig,ax,figpad)
% figureformat(fig,ax,figpad)
% Formats input figure to reduce white space padding around plot to some
% value set by figpad (0.005 typical).  Preserves plot's aspect ratio and
% resizes to have a 500 pixel width.

outerpos = ax.OuterPosition;
ti = [0 0 0 0];
left        = outerpos(1) + ti(1) + figpad;
bottom      = outerpos(2) + ti(2) + figpad;
ax_width    = outerpos(3) - ti(1) - ti(3) - 2*figpad;
ax_height   = outerpos(4) - ti(2) - ti(4) - 2*figpad;

set(ax,'Visible','on')
set(ax,'XAxisLocation','origin')
set(ax,'YAxisLocation','origin')
set(ax,'PlotBoxAspectRatio',[1 1 1])
set(ax,'TickLabelInterpreter','latex')
ax.Position = [left bottom ax_width ax_height];

ax.Position = [left bottom ax_width ax_height];

aspect = ax_width/ax_height;
fig.Position = [0 0 500 500/aspect];
box on
end
