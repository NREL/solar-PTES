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

% Obtain default or manual dimensions
if nargin == 3 % default
    width  = 14.8167;
    height = 11.1125;
elseif nargin == 4 % take manually given values
    width  = dimensions(1);
    height = dimensions(2);
else
    error('Number of inputs not implemented')
end

% Grap figure and set dimensions
figure(num)
fig = gcf;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 width height]; %default size

% Make sure that 'formats' is a cell array
if ~iscell(formats)
    formats = {formats};
end

% Export figure using saveas
for i=1:length(formats)
    % Special commands are needed for .emf files
    if strcmp(formats{i},'emf')
        switch computer
            case 'GLNXA64'
                % Save first as a temporary .svg file
                temp_name = [file_name,'_temp']; 
                saveas(fig,temp_name,'svg')
                % Use inkscape to convert .svg file to .emf
                command = ['inkscape --file ',temp_name,'.svg --export-emf ',file_name,'.emf'];                
                [~,~] = system(command);
                % Remove temporary .svg file
                [~,~] = system(['rm ',temp_name,'.svg']);
                
            case 'PCWIN64'
                warning(strcat('Verify that exporting figures to .emf',...
                    'works properly in the Windows version of Matlab'));
                saveas(fig,file_name,formats{i})
        end
        
    else
        saveas(fig,file_name,formats{i})        
    end    
end

end

