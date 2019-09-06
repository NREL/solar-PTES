function save_fig( num, file_name, width, height, mode )
% Set figure size and save in .fig and .eps formats
% file_name is the file path, e.g. './results/my_figure'
% Usage example: save_fig(1,'./results/my_figure',0,0,0);

if mode == 0 % default
    width  = 14.8167;
    height = 11.1125;
elseif mode ==1 % take manually given values
else
    error('mode not implemented')
end

figure(num)
fig = gcf;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 width height]; %default size

saveas(fig,file_name,'fig')
%saveas(fig,file_name,'epsc')
saveas(fig,file_name,'svg')

end

