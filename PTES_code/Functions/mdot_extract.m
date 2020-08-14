function mdot_multirun = mdot_extract(Npnt,Ncrv,W0)

% Go and look at the power output from a previous multi-run. Then
% recalculate what the mass flow rate should be in each case to give the
% desired mass flow rate.

% W0 is the desired power output during discharge, W

% Create a matrix by extracting and arranging the values of the variable
% 'var_name' from the different files inside the Multi_run folder.
mdot_multirun   = zeros(Npnt,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);
        S1 = load(filename,'E_net_dis');
        S2 = load(filename,'t_dis');
        S3 = load(filename,'mdot');

        Wact = S1.E_net_dis / S2.t_dis ; % Actual average power output, W
        fac  = W0 / Wact ;
        
        if fac < 0
            fac = 1 ;
        end

        mdot_multirun(ipnt,icrv)  = S3.mdot * fac;
    end
end

end