classdef LIB_THERM
    %LIB_THERM Library of thermophysical properties
    %   Each object of this class corresponds to a specific fluid.
    
    properties
        name
        list
        
        pcrit
        
        X
        Y
        F
        TAB
        
        XL
        YL
        TABL
        XLr
        YLr
        FLr
        
        XG
        YG
        TABG
        XGr
        YGr
        FGr
        
        XT
        YT
        TABT
        XTr
        YTr
        FTr
        
        xL
        yL
        tabL
        fL
        
        xG
        yG
        tabG
        fG
        
        dxL
        dxG
    end
    
    methods
        
        function obj = LIB_THERM(name,inputs,num)
            % Create a table of thermophysical properties. Two
            % interpolation tables are created: One regular table that
            % covers the whole space, and a smaller, denser table, which is
            % positioned around the two-phase region and superimposed on
            % the first one. This second grid is divided in three
            % sub-grids, one on the liquid side of the saturation curve
            % (L), one on the gas side (G), and one on top of the critical
            % pressure line (T).
            
            % Create fluid with CoolProp and obtain fluid handle
            fluid  = fluid_class(name,'WF','CP','HEOS',1,1);
            handle = fluid.handle;
            
            switch name
                case 'Water'
                    % Set minimum and maximum temperatures, pressures and enthalpies
                    Tmin  = CP1(0,0,0,'Tmin',fluid.handle);
                    Tmax  = 1500;
                    pmin  = CP1(0,0,0,'pmin',fluid.handle);
                    pmax  = 1000e5;
                    pcrit = CP1(0,0,0,'Pcrit',handle);
                    hmin  = CP1('PT_INPUTS',pmax,Tmin,'H',fluid.handle);
                    hmax  = CP1('PT_INPUTS',pmin,Tmax,'H',fluid.handle);
                otherwise
                    error('not implemented')
            end
            
            switch inputs
                case 'HmassP_INPUTS'
                    % Set list of properties to extract
                    list = {'T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'};
                    
                    % Set arrays for base grid
                    y_array = logspace(log10(pmin),log10(pmax),num);
                    x_array = linspace(hmin,hmax,num);
                    
                    % Set number of points for sub-grids around saturation curve
                    nsy = round(2*num);
                    nsx = round(nsy/10);
                    
                    % Obtain properties along the saturation curve. Use an array of
                    % pressure points between Pmin and Pcrit. 'L' denotes the
                    % liquid side of the curve (vapour quality = 0.0). 'G' denotes
                    % the gas side (vapour quality = 1.0).
                    ns    = 1e3;
                    
                    fL=@(P) CP1('PQ_INPUTS',P,0.0*ones(size(P)),'H',handle);
                    fG=@(P) CP1('PQ_INPUTS',P,1.0*ones(size(P)),'H',handle);
                    [yL,xL] = discrete_curve(fL,pmin,pcrit,pcrit*1.001,ns,1);
                    [yG,xG] = discrete_curve(fG,pcrit,pmin,pcrit*1.001,ns,0);
                    
                    % Create a regular, rectangular grid over the whole domain
                    [X,Y] = meshgrid(x_array,y_array);
                    
                    % Create enthalpy and pressure arrays for two curves that
                    % surround the saturation dome and define the Middle region.
                    % 'dH' is the separation between the saturation curve and the
                    % outer curves on the x axis. 'dP' is the separation on the y
                    % axis. one point at Pcrit.
                    dxL  = (xG(1) - xL(1))*0.15;
                    dxG  = dxL*0.5;
                    
                    % Use the outer curves to create a grid for the Middle region.
                    % YM is a matrix containing the y coordinates of each grid
                    % point. It is created using the yM sub-array which goes from
                    % Pmin to Pcrit/2. XM is a matrix containing the x coordinates
                    % of each grid point. It has two branches. The left-side
                    % branch, contained between xL (x coordinates of liquid
                    % saturation curve) and xL - dxL. And the right-side branch,
                    % contained between xG (x coordinates of gas saturation curve)
                    % and xG + dHx.
                    yM  = logspace(log10(pmin),log10(pcrit),nsy)';
                    xML = interp1(yL,xL,yM,'spline');
                    xMG = interp1(yG,xG,yM,'spline');
                    YL = yM.*ones([length(yM),nsx]);
                    YG = yM.*ones([length(yM),ceil(nsx*dxG/dxL)]);
                    XL = zeros(size(YL));
                    XG = zeros(size(YG));
                    for iM=1:length(yM)
                        XL(iM,:) = linspace(xML(iM)-dxL,xML(iM),size(XL,2));
                        XG(iM,:) = linspace(xMG(iM),xMG(iM)+dxG,size(XG,2));
                    end
                    
                    % Top part
                    dP  = pcrit/2;
                    nT  = round(nsy/2);
                    xTa = CP1('PQ_INPUTS',pcrit-dP,0.0,'H',handle);
                    xTb = CP1('PQ_INPUTS',pcrit-dP,1.0,'H',handle);
                    xT  = linspace(xTa,xTb,nT);
                    yT  = zeros(size(xT));
                    XT  = xT.*ones([nsx*2,nT]);
                    YT  = zeros(size(XT));
                    for i=1:nT
                        if xT(i) <= xL(end)
                            fun = @(P) CP1('PQ_INPUTS',P,0.0*ones(size(P)),'H',handle) - xT(i);
                        else
                            fun = @(P) CP1('PQ_INPUTS',P,1.0*ones(size(P)),'H',handle) - xT(i);
                        end
                        yT(i)   = fzero(fun,[pcrit/3,pcrit]);
                        YT(:,i) = logspace(log10(yT(i)),log10(pcrit+dP),nsx*2)';
                    end
                    
                    %{
                    figure(2)
                    semilogy(xL,yL,xG,yG,xL-dxL,yL,xG+dxG,yG); hold on;
                    semilogy(XL,YL,'k.');
                    semilogy(XG,YG,'k.');
                    semilogy(XT,YT,'k.');
                    hold off;
                    xlabel('Enthalpy, J/kg/K')
                    ylabel('Pressure, bar')
                    keyboard
                    %}
                    
                    % Obtain thermophysical properties at each grid point.
                    np    = length(list);
                    TAB   = zeros([size(X),np]);
                    TABL  = zeros([size(XL),np]);
                    TABG  = zeros([size(XG),np]);
                    TABT  = zeros([size(XT),np]);
                    tabL  = zeros([size(xL),np]);
                    tabG  = zeros([size(xG),np]);
                    for ip=1:np
                        prop = list{ip};
                        
                        switch prop
                            case {'H','T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'}
                                
                                for ic=1:size(X,2) %fill column by column
                                    TAB(:,ic,ip) = CP1(inputs,X(:,ic),Y(:,ic),prop,handle);
                                end
                                
                                for ic=1:size(XL,2)
                                    TABL(:,ic,ip) = CP1(inputs,XL(:,ic),YL(:,ic),prop,handle);
                                end
                                
                                for ic=1:size(XG,2)
                                    TABG(:,ic,ip) = CP1(inputs,XG(:,ic),YG(:,ic),prop,handle);
                                end
                                
                                for ic=1:size(XT,2)
                                    TABT(:,ic,ip) = CP1(inputs,XT(:,ic),YT(:,ic),prop,handle);
                                end
                                
                                tabL(:,ip) = CP1('PQ_INPUTS',yL,0.0*ones(size(yL)),prop,handle);
                                tabG(:,ip) = CP1('PQ_INPUTS',yG,1.0*ones(size(yG)),prop,handle);
                                
                            otherwise
                                error('not implemented')
                        end
                    end
                    
                    
                    % Save parameters inside object
                    obj.name  = name;
                    obj.list  = list;
                    
                    obj.TAB   = TAB;
                    obj.X     = X;
                    obj.Y     = Y;
                    
                    obj.pcrit = pcrit;
                    
                    obj.TABL  = TABL;
                    obj.XL    = XL;
                    obj.YL    = YL;
                    
                    obj.TABG  = TABG;
                    obj.XG    = XG;
                    obj.YG    = YG;
                    
                    obj.TABT  = TABT;
                    obj.XT    = XT;
                    obj.YT    = YT;
                    
                    obj.tabL  = tabL;
                    obj.xL    = xL;
                    obj.yL    = yL;
                    
                    obj.tabG  = tabG;
                    obj.xG    = xG;
                    obj.yG    = yG;
                    
                    obj.dxL = dxL;
                    obj.dxG = dxG;
                    
                case 'PQ_INPUTS'
                    
                    % Set list of properties to extract
                    list = {'H','T','S','D','U','C','CONDUCTIVITY','VISCOSITY','PRANDTL'};
                    
                    % Obtain properties along the saturation curve. Use an array of
                    % pressure points between Pmin and Pcrit.
                    xs = logspace(log10(pmin),log10(pcrit),num)';
                    
                    % Obtain thermophysical properties at each grid point.
                    np    = length(list);
                    tabL  = zeros([size(xs),np]);
                    tabG  = zeros([size(xs),np]);
                    for ip=1:np
                        prop = list{ip};
                        
                        switch prop
                            case {'H','T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'}
                                
                                tabL(:,ip) = CP1('PQ_INPUTS',xs,0.0*ones(size(xs)),prop,handle);
                                tabG(:,ip) = CP1('PQ_INPUTS',xs,1.0*ones(size(xs)),prop,handle);
                                
                            otherwise
                                error('not implemented')
                        end
                    end
                    
                    % Save parameters inside object
                    obj.name  = name;
                    obj.list  = list;
                    
                    obj.tabL  = tabL;
                    obj.xL    = xs;
                    
                    obj.tabG  = tabG;
                    obj.xG    = xs;
                    
                otherwise
                    error('not implemented')
            end
        end
        
        function [obj, results] = RLIB(obj, inputs, Xq, Yq, out)
            %RLIB Read library of thermophysical properties.
            
            switch inputs
                case 'HP'
                    
                    % Create gridded interpolants
                    if isempty(obj.F)
                        np     = length(obj.list);
                        obj.F  = cell([1 np]);
                        obj.fL = cell([1 np]);
                        obj.fG = cell([1 np]);
                        meth   = 'linear';
                        for ip=1:np
                            prop = obj.list{ip};
                            switch prop
                                case {'H','T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'}
                                    vL = obj.tabL(:,ip);
                                    vG = obj.tabG(:,ip);
                                    obj.fL{ip} = griddedInterpolant(obj.yL,vL,meth);
                                    obj.fG{ip} = griddedInterpolant(obj.yG,vG,meth);
                                    obj.F{ip}  = griddedInterpolant(obj.X',obj.Y',obj.TAB(:,:,ip)',meth);
                                    
                                otherwise
                                    error('not implemented')
                            end
                        end
                        %{
                method_gi = 'linear';
                [obj.XLr,obj.YLr,obj.FLr] = normalise_grid(obj.XL,obj.YL,obj.TABL,method_gi);
                [obj.XGr,obj.YGr,obj.FGr] = normalise_grid(obj.XG,obj.YG,obj.TABG,method_gi);
                [obj.XTr,obj.YTr,obj.FTr] = normalise_grid(obj.XT,obj.YT,obj.TABT,method_gi);
                        %}
                    end
                    
                    % Obtain logical array corresponding to indices of the
                    % properties 'out' inside obj.list
                    ind = zeros(size(obj.list));
                    for io = 1:length(out)
                        cond = strcmp(out{io},obj.list);
                        if ~any(cond)
                            error('invalid property selected')
                        elseif any(ind & cond)
                            error('repeated property')
                        else
                            ind = ind | cond;
                        end
                    end
                    
                    % Determine which query points fall within each region
                    method_sat = 'spline';
                    XLq = interp1(obj.yL,obj.xL,Yq,method_sat);
                    XGq = interp1(obj.yG,obj.xG,Yq,method_sat);
                    SCRITq = Yq >= obj.pcrit;        % Super-critical
                    SCOOLq = Xq <= XLq & ~SCRITq;    % Sub-cooled
                    SHEATq = Xq >= XGq & ~SCRITq;    % Super-heated
                    OUTq   = SCRITq | SCOOLq | SHEATq; % Outside 2-phase region
                    INq    = ~OUTq;                    % Inside  2-phase region
                    
                    % Determine query points within intermediate regions
                    INLq = XLq - obj.dxL <= Xq & SCOOLq & Yq <= obj.YL(end,1);
                    INGq = Xq <= XGq + obj.dxG & SHEATq & Yq <= obj.YG(end,1);
                    INTq = Xq >= obj.XT(1,1) & Xq <= obj.XT(1,end) & Yq >= obj.YT(1,1) & Yq <= obj.YT(end,1) & OUTq;
                    
                    %{
            figure(3)
            semilogy(obj.X,obj.Y,'k.',obj.XL,obj.YL,'k.',obj.XG,obj.YG,'k.'); hold on;
            semilogy(Xq(INLq),Yq(INLq),'rs')
            semilogy(Xq(INGq),Yq(INGq),'bs')
            semilogy(Xq(INTq),Yq(INTq),'ys')
            hold off;
            keyboard
                    %}
                    
                    VsatLq = zeros([size(XLq),length(out)]);
                    VsatGq = zeros([size(XGq),length(out)]);
                    for io = 1:length(out)
                        ind = strcmp(out{io},obj.list);
                        VsatLq(:,:,ind) = obj.fL{ind}(Yq);
                        VsatGq(:,:,ind) = obj.fG{ind}(Yq);
                    end
                    
                    % Fill in table of vapour quality
                    Q      = ones(size(Xq)).*(-1);
                    Q(INq) = (Xq(INq) - XLq(INq))./(XGq(INq) - XLq(INq));
                    
                    %{
            % Obtain query points for normalised (rectangular) grids
            [xqLr,yqLr] = map_grid(obj.XL,obj.YL,Xq(INLq),Yq(INLq),'h');
            [xqGr,yqGr] = map_grid(obj.XG,obj.YG,Xq(INGq),Yq(INGq),'h');
            [xqTr,yqTr] = map_grid(obj.XT,obj.YT,Xq(INTq),Yq(INTq),'v');
                    %}
                    
                    % Interpolate points outside saturation dome
                    results = zeros([size(Xq),length(out)]);
                    for io = 1:length(out)
                        ind = strcmp(out{io},obj.list);
                        Vq  = zeros(size(Xq));
                        
                        Vq(OUTq) = obj.F{ind}(Xq(OUTq),Yq(OUTq));
                        
                        %{
                Vq(INLq) = obj.FLr{ind}(xqLr,yqLr);
                Vq(INGq) = obj.FGr{ind}(xqGr,yqGr);
                %Vq(INTq) = obj.FTr{ind}(xqTr,yqTr);
                        %}
                        
                        %%{
                        VL = obj.TABL(:,:,ind);
                        Vq(INLq) = rtab_nest(obj.XL,obj.YL,VL,Xq(INLq),Yq(INLq),'h');
                        
                        VG = obj.TABG(:,:,ind);
                        Vq(INGq) = rtab_nest(obj.XG,obj.YG,VG,Xq(INGq),Yq(INGq),'h');
                        
                        VT = obj.TABT(:,:,ind);
                        Vq(INTq) = rtab_nest(obj.XT,obj.YT,VT,Xq(INTq),Yq(INTq),'v');
                        %}
                        
                        
                        % Compute points inside saturation dome using correlation
                        % with vapour quality
                        VsatL = VsatLq(:,:,ind);
                        VsatG = VsatGq(:,:,ind);
                        switch out{io}
                            case 'Q'
                                Vq = Q;
                            case 'T'
                                Vq(INq) = VsatL(INq);
                            case {'S','U'}
                                Vq(INq) = Q(INq).*VsatG(INq) + (1-Q(INq)).*VsatL(INq);
                            case 'D'
                                Vq(INq) = 1./(Q(INq)./VsatG(INq) + (1-Q(INq))./VsatL(INq));
                            case {'C','CONDUCTIVITY','VISCOSITY','PRANDTL'}
                                Vq(INq) = NaN(size(Q(INq)));
                            otherwise
                                error('not implemented');
                        end
                        
                        results(:,:,io) = Vq;
                    end
                    
                otherwise
                    error('not implemented')
            end
            
            function Vq = rtab_nest(X,Y,V,xq,yq,mode)
                
                switch mode
                    case 'h'
                        % Horizontally adapted grid (constant vertical axis)
                        
                        % Select a 1D array from the Y matrix which is
                        % representative of the selected region. Store as a
                        % 1-column array.
                        Y1D  = Y(:,1);
                        ymin = min(Y1D);
                        yr   = Y1D(2)/Y1D(1);
                        numy = length(Y1D);
                        
                        % Obtain the y-indices (and values) around the yq query
                        % points.
                        iy1 = floor(1 + log(yq/ymin)./log(yr));
                        iy1(iy1<=1) = 1;
                        iy1(iy1>=numy-1) = numy - 1;
                        iy2 = iy1 + 1;
                        y1  = ymin*yr.^(iy1 - 1);
                        y2  = ymin*yr.^(iy2 - 1);
                        
                        % Obtain the x-indices (and values) around the xq query
                        % points. In this case, the values of xmin and dx vary for
                        % each point.
                        % x-indices for iy1 (a)
                        xmin_a = min(X(iy1,:),[],2);
                        dx_a   = X(iy1,2)-X(iy1,1);
                        numx = length(X(1,:));
                        ix1_a = floor(1 + (xq-xmin_a)./dx_a);
                        ix1_a(ix1_a<=1) = 1;
                        ix1_a(ix1_a>=numx) = numx - 1;
                        ix2_a = ix1_a + 1;
                        x1_a  = xmin_a + (ix1_a - 1).*dx_a;
                        x2_a  = xmin_a + (ix2_a - 1).*dx_a;
                        % x-indices for iy2 (b)
                        xmin_b = min(X(iy2,:),[],2);
                        dx_b   = X(iy2,2)-X(iy2,1);
                        ix1_b = floor(1 + (xq-xmin_b)./dx_b);
                        ix1_b(ix1_b<=1) = 1;
                        ix1_b(ix1_b>=numx) = numx - 1;
                        ix2_b = ix1_b + 1;
                        x1_b  = xmin_b + (ix1_b - 1).*dx_b;
                        x2_b  = xmin_b + (ix2_b - 1).*dx_b;
                        
                        % Extract coefficients for bilinear interpolation formula
                        % (i.e. values of V on corners around query point).
                        Q11  = V(sub2ind(size(V),iy1,ix1_a));
                        Q12  = V(sub2ind(size(V),iy2,ix1_b));
                        Q21  = V(sub2ind(size(V),iy1,ix2_a));
                        Q22  = V(sub2ind(size(V),iy2,ix2_b));
                        
                        % Apply bilinear interpolation formula
                        Vy1 = ((x2_a-xq).*Q11 + (xq-x1_a).*Q21)./(x2_a-x1_a);
                        Vy2 = ((x2_b-xq).*Q12 + (xq-x1_b).*Q22)./(x2_b-x1_b);
                        Vq  = ((y2-yq).*Vy1 + (yq-y1).*Vy2)./(y2-y1);
                        
                    case 'h2'
                        % Horizontally adapted grid (constant vertical axis)
                        % use spline method!
                        
                        % Select a 1D array from the Y matrix which is
                        % representative of the selected region. Store as a
                        % 1-column array.
                        Y1D  = Y(:,1);
                        ymin = min(Y1D);
                        yr   = Y1D(2)/Y1D(1);
                        numy = length(Y1D);
                        
                        % Obtain the y-indices (and values) around the yq query
                        % points.
                        iy1 = floor(1 + log(yq/ymin)./log(yr));
                        iy1(iy1<=2) = 2;
                        iy1(iy1>=numy-2) = numy - 2;
                        iy0 = iy1 - 1;
                        iy2 = iy1 + 1;
                        iy3 = iy1 + 2;
                        
                        method = 'spline';
                        Vq  = zeros(size(xq));
                        for iq=1:length(xq)
                            Vy0 = interp1(X(iy0(iq),:),V(iy0(iq),:),xq(iq),method);
                            Vy1 = interp1(X(iy1(iq),:),V(iy1(iq),:),xq(iq),method);
                            Vy2 = interp1(X(iy2(iq),:),V(iy2(iq),:),xq(iq),method);
                            Vy3 = interp1(X(iy3(iq),:),V(iy3(iq),:),xq(iq),method);
                            Vq(iq) = interp1(Y1D(iy0(iq):iy3(iq)),[Vy0;Vy1;Vy2;Vy3],yq(iq),method);
                        end
                        
                    case 'v'
                        % Vertically adapted grid (constant horizontal axis)
                        
                        % Select a 1D array from the X matrix which is
                        % representative of the selected region. Store as a
                        % 1-column array.
                        X1D   = X(1,:)';
                        xmin  = min(X1D);
                        dx    = X1D(2) - X1D(1);
                        numx  = length(X1D);
                        
                        % Obtain the x-indices (and values) around the Xq query
                        % points
                        ix1 = floor(1 + (xq-xmin)./dx);
                        ix1(ix1<=1) = 1;
                        ix1(ix1>=numx) = numx - 1;
                        ix2 = ix1 + 1;
                        x1  = X1D(ix1);%xmin + (ix1 - 1).*dx
                        x2  = X1D(ix2);%xmin + (ix2 - 1).*dx
                        
                        % Obtain the y-indices (and values) around the yq query
                        % points. In this case, the values of ymin and yr vary for
                        % each point.
                        % y-indices for ix1 (a)
                        ymin_a = min(Y(:,ix1))';
                        yr_a   = (Y(2,ix1)./Y(1,ix1))';
                        numy   = length(Y(:,1));
                        iy1_a = floor(1 + log(yq./ymin_a)./log(yr_a));
                        iy1_a(iy1_a<=1) = 1;
                        iy1_a(iy1_a>=numy) = numy - 1;
                        iy2_a = iy1_a + 1;
                        y1_a  = Y(sub2ind(size(Y),iy1_a,ix1));%ymin_a.*yr_a.^(iy1_a - 1);
                        y2_a  = Y(sub2ind(size(Y),iy2_a,ix1));%ymin_a.*yr_a.^(iy2_a - 1);
                        % x-indices for ix2 (b)
                        ymin_b = min(Y(:,ix2))';
                        yr_b   = (Y(2,ix2)./Y(1,ix2))';
                        iy1_b = floor(1 + log(yq./ymin_b)./log(yr_b));
                        iy1_b(iy1_b<=1) = 1;
                        iy1_b(iy1_b>=numy) = numy - 1;
                        iy2_b = iy1_b + 1;
                        y1_b  = Y(sub2ind(size(Y),iy1_b,ix2));%ymin_b.*yr_b.^(iy1_b - 1);
                        y2_b  = Y(sub2ind(size(Y),iy2_b,ix2));%ymin_b.*yr_b.^(iy2_b - 1);
                        
                        % Extract coefficients for bilinear interpolation formula
                        % (i.e. values of V on corners around query point).
                        Q11  = V(sub2ind(size(V),iy1_a,ix1));
                        Q12  = V(sub2ind(size(V),iy2_a,ix1));
                        Q21  = V(sub2ind(size(V),iy1_b,ix2));
                        Q22  = V(sub2ind(size(V),iy2_b,ix2));
                        
                        % Apply bilinear interpolation formula
                        Vx1 = ((y2_a-yq).*Q11 + (yq-y1_a).*Q12)./(y2_a-y1_a);
                        Vx2 = ((y2_b-yq).*Q21 + (yq-y1_b).*Q22)./(y2_b-y1_b);
                        Vq  = ((x2-xq).*Vx1 + (xq-x1).*Vx2)./(x2-x1);
                        
                    case 'v2'
                        % Vertically adapted grid (constant horizontal axis)
                        
                        % Select a 1D array from the X matrix which is
                        % representative of the selected region. Store as a
                        % 1-column array.
                        X1D   = X(1,:)';
                        xmin  = min(X1D);
                        dx    = X1D(2) - X1D(1);
                        numx  = length(X1D);
                        
                        % Obtain the x-indices (and values) around the Xq query
                        % points
                        ix1 = floor(1 + (xq-xmin)./dx);
                        ix1(ix1<=2) = 2;
                        ix1(ix1>=numx-2) = numx - 2;
                        ix0 = ix1 - 1;
                        ix2 = ix1 + 1;
                        ix3 = ix1 + 2;
                        
                        method = 'spline';
                        Vq  = zeros(size(yq));
                        for iq=1:length(yq)
                            Vx0 = interp1(Y(:,ix0(iq)),V(:,ix0(iq)),yq(iq),method);
                            Vx1 = interp1(Y(:,ix1(iq)),V(:,ix1(iq)),yq(iq),method);
                            Vx2 = interp1(Y(:,ix2(iq)),V(:,ix2(iq)),yq(iq),method);
                            Vx3 = interp1(Y(:,ix3(iq)),V(:,ix3(iq)),yq(iq),method);
                            Vq(iq) = interp1(X1D(ix0(iq):ix3(iq)),[Vx0;Vx1;Vx2;Vx3],xq(iq),method);
                        end
                        
                    case 'g'
                        % General, rectangular grid
                        
                        % Select a 1D array from the Y matrix. Store as a
                        % 1-column array.
                        Y1D   = Y(:,1);
                        ymin  = min(Y1D);
                        yr    = Y1D(2)/Y1D(1);
                        numy  = length(Y1D);
                        
                        % Obtain the y-indices (and values) around the Yq query
                        % points.
                        iy1 = floor(1 + log(yq/ymin)./log(yr));
                        iy1(iy1<=1) = 1;
                        iy1(iy1>=numy) = numy - 1;
                        iy2 = iy1 + 1;
                        y1  = ymin*yr.^(iy1 - 1);
                        y2  = ymin*yr.^(iy2 - 1);
                        
                        % Select a 1D array from the X matrix. Store as a
                        % 1-column array.
                        X1D   = X(1,:)';
                        xmin  = min(X1D);
                        dx    = X1D(2) - X1D(1);
                        numx  = length(X1D);
                        
                        % Obtain the x-indices (and values) around the Xq query
                        % points
                        ix1 = floor(1 + (xq-xmin)./dx);
                        ix1(ix1<=1) = 1;
                        ix1(ix1>=numx) = numx - 1;
                        ix2 = ix1 + 1;
                        x1  = xmin + (ix1 - 1).*dx;
                        x2  = xmin + (ix2 - 1).*dx;
                        
                        % Extract coefficients for bilinear interpolation formula
                        % (i.e. values of V on corners around query point).
                        Q11  = V(sub2ind(size(V),iy1,ix1));
                        Q12  = V(sub2ind(size(V),iy2,ix1));
                        Q21  = V(sub2ind(size(V),iy1,ix2));
                        Q22  = V(sub2ind(size(V),iy2,ix2));
                        
                        % Apply bilinear interpolation formula
                        Vy1 = ((x2-xq).*Q11 + (xq-x1).*Q21)./(x2-x1);
                        Vy2 = ((x2-xq).*Q12 + (xq-x1).*Q22)./(x2-x1);
                        Vq  = ((y2-yq).*Vy1 + (yq-y1).*Vy2)./(y2-y1);
                        
                    otherwise
                        error('not implemented')
                end
            end
        end
        
        function [Xr,Yr,Fr] = normalise_grid(X,Y,TAB,method)
            % Maps a semi-rectangular grid into a regular grid
            % contained between [0,1], and creates gridded
            % interpolants.
            xr = linspace(0,1,size(X,2));
            yr = linspace(0,1,size(Y,1));
            [Xr,Yr] = meshgrid(xr,yr);
            
            nV  = size(TAB,3);
            Fr  = cell([1 nV]);
            for iV=1:nV
                V     = TAB(:,:,iV);
                Fr{iV} = griddedInterpolant(Xr',Yr',V',method);
            end
            
        end
        
        function [xqr,yqr] = map_grid(X,Y,xq,yq,mode)
            % Mapping principle:
            % Xr  = (X - X(:,1))./(X(:,end) - X(:,1));
            % Yr  = (Y - Y(1,:))./(Y(end,:) - Y(1,:));
            
            %Y  = log10(Y);
            %yq = log10(yq);
            %y1 = log10(y1);
            
            method = 'linear';
            
            switch mode
                case 'h'
                    % Horizontally adapted grid (constant vertical
                    % axis)
                    
                    % Select a 1D array from the Y matrix. Store as a
                    % 1-column array.
                    Y1D = Y(:,1);
                    
                    % Obtain initial and final positions of the x-axis
                    % grid corresponding to each yq point
                    xa  = interp1(Y1D,X(:,1),yq,method);
                    xb  = interp1(Y1D,X(:,end),yq,method);
                    
                    figure(4)
                    semilogy(X,Y,'k.'); hold on;
                    semilogy(xq,yq,'r.');
                    semilogy(xa,yq,'b.');
                    semilogy(xb,yq,'y.');
                    hold off;
                    
                    % Map xq and yq query points onto normalised
                    % (regular) grid
                    xqr = (xq - xa)./(xb - xa);
                    yqr = (yq - Y1D(1))./(Y1D(end) - Y1D(1));
                    
                    keyboard
                    
                case 'v'
                    % Vertically adapted grid (constant horizontal
                    % axis)
                    
                    % Select a 1D array from the X matrix. Store as a
                    % 1-column array.
                    X1D = X(1,:);
                    
                    % Obtain initial and final positions of the y-axis
                    % grid corresponding to each xq point
                    ya  = interp1(X1D,Y(1,:),xq,method);
                    yb  = interp1(X1D,Y(end,:),xq,method);
                    
                    % Map xq and yq query points onto normalised
                    % (regular) grid
                    xqr = (xq - X1D(1))./(X1D(end) - X1D(1));
                    yqr = (yq - ya)./(yb - ya);
                    
                otherwise
                    error('not implemented')
            end
            
        end
    end
end

