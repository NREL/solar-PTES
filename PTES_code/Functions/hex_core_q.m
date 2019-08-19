function [QT] = hex_core_q(TvH, CpvH, hvH, TvC, CpvC, hvC, mH, mC, TC1, TH2, eff)
% Obtain QMAX by finding temperature distribution for which DTp=0

mode = 0;
switch mode
    case 0 % Use fact that mC*CpC = mH*CpH at pinch point (faster)
        [ pp, np ] = find_pinch_points(TvH,CpvH,TvC,CpvC,mH,mC);
        
        % Possible values of Qp; total heat transfer for defined pinch point temperature(s)
        Qp = zeros(1,np);
        for i0=1:np % np is the number of possible pinch points within the HEX
            Tp  = pp(i0); % pp is an array containing the possible pinch point temperature(s)
            QC = mC*(rtab1(TvC,hvC,Tp,0)-rtab1(TvC,hvC,TC1,0)); % mC*DhC from TC1 to Tp
            QH = mH*(rtab1(TvH,hvH,TH2,0)-rtab1(TvH,hvH,Tp,0)); % mH*DhH from Tp to TH2
            Qp(i0) = QC + QH;
        end
        
        % Compute total maximum heat transfer for the separate streams, QCT and QHT
        QCT = mC*(rtab1(TvC,hvC,TH2,0) - rtab1(TvC,hvC,TC1,0)) ;
        QHT = mH*(rtab1(TvH,hvH,TH2,0) - rtab1(TvH,hvH,max([TC1,min(TvH)]),0)); %account for minimum possible TH
        
        % QMAX has to be the minimum between Qp, QC and QH:
        QMAX = min([Qp,QCT,QHT]);
        

    case 1 % Compute temperature distributions and DT at each point (slower)
        
        % Compute TH1 for which DT_min = 0
        n = length(TvH);
        f = @(TH1) min(abs(compute_DT(TH1,TH2,TvH,hvH,TC1,TvC,hvC,n)));
        [TH1,~,~,~,~] = golden_search(f,TC1,TH2,0.0001,'Min',100);
        
        % Compute QMAX for the TH1 value found
        hH1 = rtab1(TvH,hvH,TH1,0);
        hH2 = rtab1(TvH,hvH,TH2,0);
        QMAX = mH*(hH2 - hH1);
        
    otherwise
        error('not implemented')
end

% Actual heat transfer is QMAX times the given effectiveness
QT = QMAX*eff;

end

function DT = compute_DT(TH1,TH2,TvH,hvH,TC1,TvC,hvC,n)

% Compute temperature distribution of hot stream
TH  = linspace(TH1,TH2,n)';
hH  = rtab1(TvH,hvH,TH,0);
hH1 = hH(1);

% Compute temperature distribution of cold stream
QS  = hH - hH1; % cummulative heat transfer
hC1 = rtab1(TvC,hvC,TC1,0);
hC  = hC1 + QS;
TC  = rtab1(hvC,TvC,hC,1);

% Compute temperature difference between the two streams
DT  = TH - TC;

end

