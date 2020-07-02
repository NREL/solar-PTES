% Set default options for plots
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',1.5)

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultAxesFontSize',14);
set(0,'defaultLegendFontSize',12);
set(0,'defaultAxesLabelFontSizeMultiplier',1.25)
set(0,'defaultAxesTitleFontSizeMultiplier',1.25)
set(0,'defaultfigurecolor','w')

% Get colors: (or manually using "uisetcolor")
c_pale_blue   = [0.3020    0.7490    0.9294];
c_dark_blue   = [0.0000    0.4510    0.7412];
c_pale_green  = [0.4706    0.6706    0.1882];
c_dark_green  = [0.1020    0.4000    0.1020];
c_pale_orange = [0.9294    0.6902    0.1294];
c_dark_orange = [0.8510    0.3294    0.1020];
c_yellow      = [1.0000    1.0000    0.0667];
c_grey        = [0.6510    0.6510    0.6510];
c_dark_grey   = [0.4000    0.4000    0.4000];

% Set arrays containing several ticks and markers for plots
sty{1} = '-'; 
sty{2} = '-.';
sty{3} = ':'; 
sty{4} = '--';
sty{5} = '-'; 
sty{6} = '-.';
sty{7} = ':'; 
sty{8} = '--';

stym{1} = 'o'; 
stym{2} = '^'; 
stym{3} = 'v'; 
stym{4} = 's'; 
stym{5} = '+'; 
stym{6} = 'x'; 
stym{7} = 'p'; 
stym{8} = '>'; 