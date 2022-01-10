function [oldAxis,newAxis] = magnifyPlot(xlm,ylm,pos,lwid,showAxisLabels)
% This function magnifies a specified rectangular area [xlm, ylm]
% in the current figure,
% and plots the magnified area on the same figure at axis defined by 
% position pos
%
% Inputs:
%       xlm = 1x2 vector specifying limits of x-axis
%       ylm = 1x2 vector specifying limits of y-axis
%       pos = 4x1 vector specifying position of new axis
%       lwid = linewidth of lines in magnified plot
%       showAxisLabels = boolean to turn axis labels on/off
%
% Example usage:
%       x = 0:10; y1 = x; y2 = 3 - 0.5*x; 
%       figure; 
%       plot(0:10,y1,'k--','linewidth',2); hold on;
%       plot(0:10,y2,'ko-','linewidth',2); hold on;
%       axis square;
%       xlm = [1.5,2.5];
%       ylm = [1.5,2.5];
%       pos = [.25 .6 .25 .25];
%
%       magnifyPlot(xlm,ylm,pos,2,false)
%
% (c) Rishabh Datta, 2021-02-12

oldAxis = gca;
ch = get(gca,'children'); % get current figure children
for ii = 1:numel(ch)
X{ii} = get(ch(ii),'xdata'); % store all data on current figure
Y{ii} = get(ch(ii),'ydata');
clr{ii} = get(ch(ii),'Color');
ls{ii} = get(ch(ii),'LineStyle');
mkr{ii} = get(ch(ii),'Marker');
mk_clr{ii} = get(ch(ii),'MarkerFaceColor');
end
% Plot magnification rectangle
v = [xlm(1) ylm(1); ...
    xlm(2) ylm(1); ...
    xlm(2) ylm(2);...
    xlm(1) ylm(2)]; % vertices of rectangle
patch('Faces',[1,2,3,4],'Vertices',v,...
    'Edgecolor','k',...
    'Facecolor','k',...
    'Facealpha',0.1,...
    'Linewidth',2,'HandleVisibility','off');% plot rectangle
axes('pos',pos) % Create axis 
for ii = 1:numel(X) % plot each child
plot(X{ii},Y{ii},'Color',clr{ii},'LineStyle',ls{ii},...
    'Marker',mkr{ii},'MarkerFaceColor',mk_clr{ii},'Linewidth',lwid); hold on;
end
xlim(xlm); ylim(ylm); % set limits of zoomed plot
set(gca,'XColor','k', 'YColor','k');
set(gca, 'LineWidth',3);
if ~showAxisLabels % turn on/off labels 
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
end
newAxis = gca;
end