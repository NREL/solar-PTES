function [ WR, WR_dis ] = work_ratio(stage,n1,n2)
%Find out the work ratio during charge and discharge

work_comp = 0;
work_exp = 0;
for i=1:n1
    if stage(1,i).w >= 0
        work_exp  = work_exp  + stage(1,i).w;
    else
        work_comp = work_comp + stage(1,i).w;
    end
end
WR = (-work_comp)/(work_exp);


work_comp = 0;
work_exp = 0;
for i=1:n2
    if stage(2,i).w >= 0
        work_exp  = work_exp  + stage(2,i).w;
    else
        work_comp = work_comp + stage(2,i).w;
    end
end
WR_dis = (work_exp)/(-work_comp);


end

