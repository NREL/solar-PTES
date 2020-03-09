function movefiles(fname)

mkdir(fname)
copyfile('./Outputs/Costs.*',fname)
copyfile('./Outputs/Losses.*',fname)
copyfile('./Outputs/T-s.*',fname)
copyfile('./Outputs/JB_RANK_INPUTS.m',fname)
copyfile('./Outputs/log.txt',fname)
copyfile('./Outputs/vars.mat',fname)

end