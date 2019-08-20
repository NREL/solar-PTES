function v = iphase_gas()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 83);
  end
  v = vInitialized;
end
