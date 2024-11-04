function [ hFig ] = plotChart( hFig, s, ucl, delay, window )
%PLOTCHART plots the monitoring statistics over time. At each time period,
%the new monitoring statistics is added to the plot.
%
% Required input arguments: 
%   hFig    : Figure handler where the new observation should be ploted. for 
%             the first observation it is set as empty.
%   s       : Monitoring statistic at the given time.
%   ucl     : Upper control limit at the given time.
%
% Optional input arguments:
%   delay   : Waithing tiem between two consecutive plots.
%   window  : Size of the window to show in the plot.
%
% I/O: [ hFig ] = plotChart( hFig, s, ucl, delay, window );
%
% The output of PLOTCHART is:
%   hFig    : Updated figure handler. Needed to do plots in loop.


if nargin<=3,
    window=20;
end
if nargin<=4,
    window=20;
end

if isempty(hFig)==1,%ishandle(hFig)==0,
    n=size(s,1);
    hFig(1)=plot([1:n]',s,'-k.');
    hold on
    hFig(2)=plot([1:n]',ucl,'-r');
    hold off
    xlabel('Sample')
    ylabel('stat')
else
    pause(delay)

    x = get(hFig,'XData');
    y = get(hFig,'YData');
    x_new=[x{1} x{1}(end)+1; x{2} x{2}(end)+1];
    y_new=[y{1} s; y{2} ucl];
    window=min(size(x_new,2),window);
    
    set(hFig(1),'XData',x_new(1,end-window+1:end),'YData',y_new(1,end-window+1:end));
    set(hFig(2),'XData',x_new(2,end-window+1:end),'YData',y_new(2,end-window+1:end));
    set(gca,'XTick',[x_new(1,end-window+1):2:x_new(end)+1])
    xlim([x_new(1,end-window+1) x_new(end)+1])
end


end

