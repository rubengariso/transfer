function [ ] = fun_savefig(h,filename,format)
% 09.10.2020
% 14.04.2021 - improvemts on figure crop and save

%for current figure use gcf as the first argument


switch format
    case {'wide'}
        base_fig_width  = 30;% cm
        base_fig_height = 11;% cm
    case {'normal'}
        base_fig_width  = 15;% cm
        base_fig_height = 11;% cm
end

fig_pos=[0 15 base_fig_width base_fig_height];%cm

% Synchronise units for screen and paper
h.Units = 'centimeters';
h.PaperUnits='centimeters';
% get all axes
allAxesInFigure = findall(h,'type','axes');
% Synchronise units of figure an axes object
for i=1:length(allAxesInFigure),
    set(allAxesInFigure(i),'Units','centimeters');
end

% set desired size
h.Position=fig_pos;

% Minimum margins around plot (as calculated by Matlab)
ti = get(allAxesInFigure(end),'TightInset');
min_margin_l = ti(1);
min_margin_b = ti(2);
min_margin_r = ti(3);
min_margin_t = ti(4);

% Set width and heigth of figure (in centimeters)
fig_width  = fig_pos(3);% h.Position(3);
fig_height = fig_pos(4);% h.Position(4);

% Add margins to the minimum margins
delta=1e-1/2;% cm
margin_l = min_margin_l + delta;
margin_b = min_margin_b + delta;
margin_r = min_margin_r + delta;
margin_t = min_margin_t + delta;

% Set dimensions of axes regions to match custom margins
left   = margin_l;
bottom = margin_b;
width  = fig_width  - margin_l - margin_r;
height = fig_height - margin_b - margin_t;

for i=1:length(allAxesInFigure),
    set(allAxesInFigure(i),'Position',[left, bottom, width, height]);
end

% Figure position on paper
% h.PaperPositionMode = 'auto';
h.PaperPositionMode = 'manual';
h.PaperPosition = [0 0 fig_width, fig_height];

% Set paper size to figure size
h.PaperSize = [fig_width, fig_height];

% set desired size again (some times it changes)
h.Position=fig_pos;

% adjust backgroud color
set(h, 'color','w');

% Save figure
print(h,[filename,'.png'],'-dpng','-r600')

% Some times the figure is not saved with the correct dimentions, 
% repeat save step
iter=1;
iter_max=10;
pp=get(h,'Position');

while sum(abs(pp-fig_pos))>0 && iter<iter_max,
    iter=iter+1;
    h.Position=fig_pos;
    print(h,[filename,'.png'],'-dpng','-r600')
    pp=get(h,'Position');
    if sum(abs(pp-fig_pos))>0 && iter>=iter_max,
        warning('saved with different position.')
    end
end
    
    
end

