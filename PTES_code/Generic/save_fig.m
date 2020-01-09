function save_fig( num, file_name, formats, dimensions )
% SAVE_FIG
%   This function saves a figure into specified formats.
%   "num" is the figure window number
%   "file_name" is the name of the file to be saved, without format.
%   "formats" is the list of formats in which to save the file in.
%   The variable "dimensions" (optional) is used to specify width and
%   length, in centimeters.
%
%   USAGE e.g.:
%   save_fig(1,'example',{'fig','epsc','svg'},[10,20])
%   save_fig(1,'example','epsc')

if nargin == 3 % default
    width  = 14.8167;
    height = 11.1125;
elseif nargin == 4 % take manually given values
    width  = dimensions(1);
    height = dimensions(2);
else
    error('Number of inputs not implemented')
end

figure(num)
fig = gcf;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 width height]; %default size

if ~iscell(formats)
    formats = {formats};
end
for i=1:length(formats)
    saveas(fig,file_name,formats{i})
end

end

