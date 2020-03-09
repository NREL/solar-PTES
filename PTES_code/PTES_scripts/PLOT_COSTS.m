figure(9)
switch Load.mode
    case {0,3,4,5,6}
        names = {'Charge machine','Discharge machine','Hot tanks','Cold tanks','Hot HX','Cold HX','Recuperators','Motor-generator'};
        Cmatrix = zeros(8,7) ;
        
        % Assign different types of costs to different places
        for ii = 1 : length(CCMP)
            Cmatrix(1,1) = Cmatrix(1,1) + CCMP(ii).cmpexp_cost.COST ;
        end
        
        for ii = 1 : length(DEXP)
            Cmatrix(2,2) = Cmatrix(2,2) + DEXP(ii).cmpexp_cost.COST ;
        end
        
        for ii = 1 : length(CEXP)
            Cmatrix(1,2) = Cmatrix(1,2) + CEXP(ii).cmpexp_cost.COST ;
        end
        
        for ii = 1 : length(DCMP)
            Cmatrix(2,1) = Cmatrix(2,1) + DCMP(ii).cmpexp_cost.COST ;
        end
        
        % Hot tank cost and hot fluid cost
        for ii = 1 : Nhot
            Cmatrix(3,3) = Cmatrix(3,3) + HT(ii).tankA_cost.COST ;    
            Cmatrix(3,4) = Cmatrix(3,4) + HT(ii).tankB_cost.COST ;    
            Cmatrix(3,5) = Cmatrix(3,5) + HT(ii).fluid_cost.COST ;    
        end
        
        % Cold tank cost and cold fluid cost
        for ii = 1 : Ncld
            Cmatrix(4,3) = Cmatrix(4,3) + CT(ii).tankA_cost.COST ;    
            Cmatrix(4,4) = Cmatrix(4,4) + CT(ii).tankB_cost.COST ;    
            Cmatrix(4,5) = Cmatrix(4,5) + CT(ii).fluid_cost.COST ;    
        end
        
        % Run through the heat exchangers and allocate them to columns
        % depending on their name
        for ii = 1 : length(HX)
           switch HX(ii).name
               case 'hot'
                   Cmatrix(5,6) = Cmatrix(5,6) + HX(ii).hx_cost.COST ;
               case 'cold'
                   Cmatrix(6,6) = Cmatrix(6,6) + HX(ii).hx_cost.COST ;
               case 'regen'
                   Cmatrix(7,6) = Cmatrix(7,6) + HX(ii).hx_cost.COST ;
               case 'rej'
                   %matrix(8,6) = matrix(8,6) + HX(ii).hx_cost.COST ;
           end
        end
        
        % OTHER
        % Motor-generator
        Cmatrix(8,7) = GEN.gen_cost.COST ;
        
        
        b = bar(Cmatrix,'stacked');
        set(gca, 'XTick', 1:8, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
        xtickangle(45)
    case 1
        error('Not implemented')
    case 2
        error('Not implemented')
end
b(1).FaceColor = c_dark_blue;
b(2).FaceColor = c_pale_blue;
b(3).FaceColor = c_pale_green;
b(4).FaceColor = c_yellow;
b(5).FaceColor = c_pale_orange;
b(6).FaceColor = c_dark_orange;
b(7).FaceColor = c_grey;
ylabel(strcat('Capital cost, \$'))
legend({'Compresors','Expanders','Source Tank','Sink Tank','Fluid','Heat exchangers','Other'},'Location','Best')


