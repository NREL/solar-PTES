function v = iphase_unknown()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 85);
  end
  v = vInitialized;
end
